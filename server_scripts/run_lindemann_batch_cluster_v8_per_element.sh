#!/bin/bash
# run_lindemann_batch_cluster_v8_per_element.sh
# 基于 v7，增加分元素Lindemann指数计算 (Pt, Sn, PtSn, PtSnO)

# ==================== 配置区 ====================

DIR_LIST_FILE="dir_list2.txt"
MAX_JOBS=16
USE_PARALLEL=true

LINDEMANN_SCRIPT="/home/scms/jychen/tools/cp2k/md/msd/lindemann_integrated_unwrap.py"
DUMP2XYZ_SCRIPT="/home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/dump2xyz_fixed.py"
SIMPLIFY_SCRIPT="/home/scms/jychen/tools/train-MLP/md-simplify2.py"
CONVERGENCE_SCRIPT="/home/scms/jychen/tools/cp2k/md/msd/analyze_convergence.py"
UNWRAP_ASE_SCRIPT="/home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/analyze_ptsn_trajectory.py"

# ✅ 新增：分元素计算脚本
PER_ELEMENT_SCRIPT="/home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/20251009/lin-new/lindemann_per_element_integrated.py"

NEED_PREPROCESS=true
DUMP2XYZ_ARGS="--extended --pbc 1 1 1"
SIMPLIFY_ARGS="1 md_ext.xyz 0.3 1.0"

USE_UNWRAP=true
UNWRAP_METHOD="ase"
UNWRAP_ELEMENTS=""

ATOM_SELECTION="name Pt or name Sn"
DO_COMPARISON=true
SKIP_EXISTING=false
CONVERGENCE_THRESHOLD=0.005

# ✅ 新增：是否启用分元素计算
ENABLE_PER_ELEMENT=true

OUTPUT_ROOT="/home/scms/jychen/tools/cp2k/md/msd/nnmD/dp-md/20251009/lin-new/collected_lindemann_cluster"

OUTPUT_MODE="clean"
SHOW_PROGRESS=true
UPDATE_INTERVAL=10

# ==================== 脚本开始 ====================

set -euo pipefail

# ==================== 读取目录列表 ====================

if [ ! -f "$DIR_LIST_FILE" ]; then
    echo "❌ 错误：目录列表文件不存在: $DIR_LIST_FILE"
    exit 1
fi

mapfile -t directories < <(
    grep -v '^[[:space:]]*$' "$DIR_LIST_FILE" | \
    sed 's/^"//; s/"$//' | \
    tr -d '\r'
)

total=${#directories[@]}

if [ "$total" -eq 0 ]; then
    echo "❌ 错误：目录列表为空"
    exit 1
fi

# ==================== 创建输出目录 ====================

mkdir -p "$OUTPUT_ROOT"

TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
RUN_ID="run_${TIMESTAMP}"

MASTER_CSV="${OUTPUT_ROOT}/lindemann_master_${RUN_ID}.csv"
COMPARISON_CSV="${OUTPUT_ROOT}/lindemann_comparison_${RUN_ID}.csv"
CONVERGENCE_CSV="${OUTPUT_ROOT}/convergence_master_${RUN_ID}.csv"
SUMMARY_LOG="${OUTPUT_ROOT}/summary_${RUN_ID}.log"
ERROR_LOG="${OUTPUT_ROOT}/error_${RUN_ID}.log"
PROGRESS_FILE="${OUTPUT_ROOT}/.progress_${RUN_ID}"

MASTER_LOCK="${OUTPUT_ROOT}/.master.lock"
COMPARISON_LOCK="${OUTPUT_ROOT}/.comparison.lock"
CONVERGENCE_LOCK="${OUTPUT_ROOT}/.convergence.lock"
PROGRESS_LOCK="${OUTPUT_ROOT}/.progress.lock"

cat > "$MASTER_CSV" << 'EOF'
目录,结构,温度(K),Lindemann指数,方法,耗时(s),时间戳
EOF

# ✅ 修改：扩展 COMPARISON_CSV 列
cat > "$COMPARISON_CSV" << 'EOF'
目录,结构,温度(K),Cluster_Unwrapped,Cluster_Wrapped,差异,差异%,Pt,Sn,PtSn,PtSnO,时间戳
EOF

cat > "$CONVERGENCE_CSV" << 'EOF'
目录,结构,温度(K),数据点,时间范围,变化率,收敛状态
EOF

echo "0" > "$PROGRESS_FILE"

# ==================== 函数定义 ====================

log_to_file() {
    local level="$1"
    shift
    local msg="[$(date '+%H:%M:%S')] [$$] $*"
    echo "$msg" >> "$SUMMARY_LOG"

    if [ "$level" = "ERROR" ]; then
        echo "$msg" >> "$ERROR_LOG"
    fi
}

log_clean() {
    local level="$1"
    local index="$2"
    shift 2
    local msg="$*"

    log_to_file "$level" "[$index/$total] $msg"

    case "$level" in
        SUCCESS)
            echo "✅ [$index/$total] $msg"
            ;;
        ERROR)
            echo "❌ [$index/$total] $msg" >&2
            ;;
        INFO)
            if [ "$OUTPUT_MODE" = "verbose" ]; then
                echo "ℹ️  [$index/$total] $msg"
            fi
            ;;
    esac
}

safe_append() {
    local lockfile="$1"
    local target="$2"
    local content="$3"

    local count=0
    while [ $count -lt 60 ]; do
        if mkdir "$lockfile" 2>/dev/null; then
            echo "$content" >> "$target"
            rmdir "$lockfile"
            return 0
        fi
        sleep 0.1
        ((count++))
    done

    echo "$content" >> "$target"
    return 1
}

update_progress() {
    local count=0
    while [ $count -lt 30 ]; do
        if mkdir "$PROGRESS_LOCK" 2>/dev/null; then
            local current=$(cat "$PROGRESS_FILE" 2>/dev/null | head -1 | tr -d '[:space:]' || echo "0")
            if [[ "$current" =~ ^[0-9]+$ ]]; then
                echo $((current + 1)) > "$PROGRESS_FILE"
            else
                echo "1" > "$PROGRESS_FILE"
            fi
            rmdir "$PROGRESS_LOCK"
            return 0
        fi
        sleep 0.1
        ((count++))
    done
    return 1
}

extract_info() {
    local dir="$1"
    local base=$(basename "$dir")
    local parent=$(basename "$(dirname "$dir")")
    local temp=$(echo "$base" | grep -oE 'T[0-9]+' | tr -d 'T')

    if [ -z "$temp" ]; then
        temp=$(echo "$parent" | grep -oE 'T[0-9]+' | tr -d 'T')
    fi

    echo "${parent}|${temp:-Unknown}"
}

check_file_valid() {
    local file="$1"
    [ -f "$file" ] && [ -s "$file" ]
}

safe_line_count() {
    local file="$1"
    if [ ! -f "$file" ]; then
        echo "0"
        return 0
    fi

    local count=$(wc -l < "$file" 2>/dev/null | head -1 | tr -d '[:space:]' || echo "0")

    if [[ "$count" =~ ^[0-9]+$ ]] && [ "$count" -ge 0 ]; then
        echo "$count"
    else
        echo "0"
    fi
}

safe_grep_count() {
    local pattern="$1"
    local file="$2"

    if [ ! -f "$file" ]; then
        echo "0"
        return 0
    fi

    local count=$(grep -c "$pattern" "$file" 2>/dev/null || echo "0")
    count=$(echo "$count" | head -1 | tr -d '[:space:]')

    if [[ "$count" =~ ^[0-9]+$ ]] && [ "$count" -ge 0 ]; then
        echo "$count"
    else
        echo "0"
    fi
}

# ==================== 单个目录处理函数 ====================

process_single_directory() {
    local dir="$1"
    local index="$2"

    local status="fail"

    exec 3>&1 4>&2
    if [ "$OUTPUT_MODE" = "clean" ]; then
        exec 1>>"$SUMMARY_LOG" 2>&1
    fi

    if [ ! -d "$dir" ]; then
        log_clean "ERROR" "$index" "目录不存在"
        echo "$status"
        return 1
    fi

    if ! cd "$dir" 2>/dev/null; then
        log_clean "ERROR" "$index" "无法进入目录"
        echo "$status"
        return 1
    fi

    IFS='|' read -r struct temp <<< "$(extract_info "$dir")"
    local basename=$(basename "$dir")

    if [ "$SKIP_EXISTING" = true ] && check_file_valid "lindemann_cluster_unwrapped.csv"; then
        if ! grep -q "$dir" "$MASTER_CSV" 2>/dev/null; then
            lindex=$(tail -n 1 lindemann_cluster_unwrapped.csv | cut -d',' -f1)
            timestamp=$(date '+%Y-%m-%d %H:%M:%S')
            safe_append "$MASTER_LOCK" "$MASTER_CSV" \
                "$dir,$struct,$temp,$lindex,skip,-,$timestamp"
        fi

        update_progress
        log_clean "INFO" "$index" "$basename (已有结果，跳过)"

        exec 1>&3 2>&4
        echo "skip"
        return 0
    fi

    log_clean "INFO" "$index" "开始: $basename ($struct, T=${temp}K)"

    LINDEMANN_INPUT=""

    if check_file_valid "sampling-simply.xyz"; then
        LINDEMANN_INPUT="sampling-simply.xyz"

    elif [ "$NEED_PREPROCESS" = true ]; then

        if check_file_valid "md_ext.xyz"; then
            if python3 "$SIMPLIFY_SCRIPT" $SIMPLIFY_ARGS &>/dev/null; then
                LINDEMANN_INPUT="sampling-simply.xyz"
            else
                log_clean "ERROR" "$index" "简化失败"
                exec 1>&3 2>&4
                echo "$status"
                return 1
            fi

        elif check_file_valid "md.xyz"; then
            if ! python "$DUMP2XYZ_SCRIPT" md.xyz md_ext.xyz $DUMP2XYZ_ARGS &>/dev/null; then
                log_clean "ERROR" "$index" "dump2xyz失败"
                exec 1>&3 2>&4
                echo "$status"
                return 1
            fi

            if python3 "$SIMPLIFY_SCRIPT" $SIMPLIFY_ARGS &>/dev/null; then
                LINDEMANN_INPUT="sampling-simply.xyz"
            else
                log_clean "ERROR" "$index" "简化失败"
                exec 1>&3 2>&4
                echo "$status"
                return 1
            fi

        else
            log_clean "ERROR" "$index" "未找到输入文件"
            exec 1>&3 2>&4
            echo "$status"
            return 1
        fi
    else
        log_clean "ERROR" "$index" "未找到sampling-simply.xyz"
        exec 1>&3 2>&4
        echo "$status"
        return 1
    fi

    if [ -z "$LINDEMANN_INPUT" ] || ! check_file_valid "$LINDEMANN_INPUT"; then
        log_clean "ERROR" "$index" "输入文件无效"
        exec 1>&3 2>&4
        echo "$status"
        return 1
    fi

    local wrapped_input="$LINDEMANN_INPUT"

    # Unwrap处理
    if [ "$USE_UNWRAP" = true ]; then
        if [ ! -f "sampling-simply_wrapped_original.xyz" ]; then
            cp "$LINDEMANN_INPUT" "sampling-simply_wrapped_original.xyz"
        fi

        case "$UNWRAP_METHOD" in
            "ase")
                unwrap_output="sampling-simply_ase_unwrapped.xyz"

                python3 "$UNWRAP_ASE_SCRIPT" "$LINDEMANN_INPUT" \
                    --surface-out "$unwrap_output" \
                    --cluster-out "test_cluster_centered.xyz" \
                    --index-out "index_zsplit.ndx" 2>&1 | tee -a unwrap_debug.log

                if [ $? -eq 0 ] && check_file_valid "$unwrap_output"; then
                    LINDEMANN_INPUT="$unwrap_output"
                    log_clean "INFO" "$index" "Unwrap成功"
                else
                    log_clean "ERROR" "$index" "Unwrap失败，使用原始文件"
                fi
                ;;
        esac
    fi

    start_time=$(date +%s)

    # 计算 Cluster (Pt+Sn) 的 Lindemann 指数
    python3 "$LINDEMANN_SCRIPT" \
        --coord "$LINDEMANN_INPUT" \
        --selection "$ATOM_SELECTION" 2>&1 | tee -a unwrap_debug.log

    calc_status=$?
    end_time=$(date +%s)
    calc_time=$((end_time - start_time))

    if [ $calc_status -ne 0 ]; then
        log_clean "ERROR" "$index" "计算失败"
        exec 1>&3 2>&4
        echo "$status"
        return 1
    fi

    if ! check_file_valid "lindemann_index_unwrapped.csv"; then
        log_clean "ERROR" "$index" "输出文件未生成"
        exec 1>&3 2>&4
        echo "$status"
        return 1
    fi

    # 重命名输出文件
    mv lindemann_index_unwrapped.csv lindemann_cluster_unwrapped.csv 2>/dev/null
    mv lindemann_per_frame_unwrapped.txt lindemann_cluster_unwrapped_per_frame.txt 2>/dev/null
    mv lindemann_integrated_unwrapped.png lindemann_cluster_unwrapped.png 2>/dev/null

    lindex_unwrap=$(tail -n 1 lindemann_cluster_unwrapped.csv | cut -d',' -f1)

    local lindex_wrapped=""
    local diff_abs=""
    local diff_pct=""

    # 对比计算：使用原始 wrapped 文件
    if [ "$DO_COMPARISON" = true ] && [ -f "sampling-simply_wrapped_original.xyz" ]; then
        python3 "$LINDEMANN_SCRIPT" \
            --coord "sampling-simply_wrapped_original.xyz" \
            --selection "$ATOM_SELECTION" 2>&1 | tee -a unwrap_debug.log

        if [ $? -eq 0 ] && check_file_valid "lindemann_index_unwrapped.csv"; then
            mv lindemann_index_unwrapped.csv lindemann_original_wrapped.csv 2>/dev/null
            mv lindemann_per_frame_unwrapped.txt lindemann_original_wrapped_per_frame.txt 2>/dev/null
            mv lindemann_integrated_unwrapped.png lindemann_original_wrapped.png 2>/dev/null

            lindex_wrapped=$(tail -n 1 lindemann_original_wrapped.csv | cut -d',' -f1)

            if [ -n "$lindex_unwrap" ] && [ -n "$lindex_wrapped" ]; then
                diff_abs=$(awk -v u="$lindex_unwrap" -v w="$lindex_wrapped" \
                    'BEGIN {printf "%.6f", u-w}')

                diff_pct=$(awk -v u="$lindex_unwrap" -v w="$lindex_wrapped" \
                    'BEGIN {
                        if(w==0) print "N/A"
                        else printf "%.2f", (u-w)/w*100
                    }')
            fi
        fi
    fi

    # ✅ 新增：分元素Lindemann指数计算
    local pt_lindex="0.000000"
    local sn_lindex="0.000000"
    local ptsn_lindex="0.000000"
    local ptsno_lindex="0.000000"

    if [ "$ENABLE_PER_ELEMENT" = true ] && check_file_valid "$PER_ELEMENT_SCRIPT"; then
        # 检查是否存在index_zsplit.ndx文件
        local index_arg=""
        if [ -f "index_zsplit.ndx" ]; then
            index_arg="--index index_zsplit.ndx"
            log_clean "INFO" "$index" "使用index_zsplit.ndx精确选择团簇原子"
        fi
        
        # 调用分元素计算脚本，输出格式：Pt,Sn,PtSn,PtSnO
        element_results=$(python3 "$PER_ELEMENT_SCRIPT" \
            --coord "$LINDEMANN_INPUT" \
            $index_arg \
            --stdout-only 2>/dev/null)

        if [ $? -eq 0 ] && [ -n "$element_results" ]; then
            IFS=',' read -r pt_lindex sn_lindex ptsn_lindex ptsno_lindex <<< "$element_results"
            log_clean "INFO" "$index" "分元素计算完成: Pt=$pt_lindex Sn=$sn_lindex"
        else
            log_clean "INFO" "$index" "分元素计算失败，使用默认值0"
        fi
    fi

    # 收敛性分析
    if check_file_valid "lindemann_cluster_unwrapped_per_frame.txt"; then
        line_count=$(safe_line_count "lindemann_cluster_unwrapped_per_frame.txt")

        if [ "$line_count" -ge 10 ]; then
            python3 "$CONVERGENCE_SCRIPT" \
                --input lindemann_cluster_unwrapped_per_frame.txt \
                --threshold $CONVERGENCE_THRESHOLD \
                --output convergence_analysis.txt 2>&1 | tee -a unwrap_debug.log

            if [ $? -eq 0 ] && check_file_valid "convergence_analysis.txt"; then
                points=$(grep "数据点数量:" convergence_analysis.txt 2>/dev/null | awk '{print $2}')
                time_start=$(grep "时间范围:" convergence_analysis.txt 2>/dev/null | awk '{print $2}')
                time_end=$(grep "时间范围:" convergence_analysis.txt 2>/dev/null | awk '{print $4}')
                time_r="${time_start}-${time_end}"
                rate=$(grep "后半段平均变化率:" convergence_analysis.txt 2>/dev/null | awk '{print $2}')
                conv_status=$(grep "结论:" convergence_analysis.txt 2>/dev/null | \
                    sed 's/.*结论://; s/✅//g; s/❌//g; s/系统//g' | xargs)

                if [ -n "$points" ] && [ -n "$conv_status" ]; then
                    safe_append "$CONVERGENCE_LOCK" "$CONVERGENCE_CSV" \
                        "$dir,$struct,$temp,$points,$time_r,$rate,$conv_status"
                fi
            fi
        fi
    fi

    # 记录结果
    if [ -n "$lindex_unwrap" ]; then
        timestamp=$(date '+%Y-%m-%d %H:%M:%S')

        safe_append "$MASTER_LOCK" "$MASTER_CSV" \
            "$dir,$struct,$temp,$lindex_unwrap,${UNWRAP_METHOD},$calc_time,$timestamp"

        # ✅ 修改：扩展 comparison CSV，增加 Pt, Sn, PtSn, PtSnO 列
        if [ "$DO_COMPARISON" = true ]; then
            safe_append "$COMPARISON_LOCK" "$COMPARISON_CSV" \
                "$dir,$struct,$temp,$lindex_unwrap,$lindex_wrapped,$diff_abs,$diff_pct,$pt_lindex,$sn_lindex,$ptsn_lindex,$ptsno_lindex,$timestamp"
        fi

        update_progress

        if [ -n "$diff_pct" ] && [ "$diff_pct" != "N/A" ]; then
            log_clean "SUCCESS" "$index" "$basename: δ(cluster)=$lindex_unwrap δ(w)=$lindex_wrapped Δ=${diff_pct}% | Pt=$pt_lindex Sn=$sn_lindex (${calc_time}s)"
        else
            log_clean "SUCCESS" "$index" "$basename: δ(cluster)=$lindex_unwrap | Pt=$pt_lindex Sn=$sn_lindex (${calc_time}s)"
        fi

        exec 1>&3 2>&4
        echo "success"
        return 0
    else
        log_clean "ERROR" "$index" "未获取到结果"
        exec 1>&3 2>&4
        echo "$status"
        return 1
    fi
}

export -f process_single_directory
export -f log_clean log_to_file
export -f safe_append update_progress extract_info check_file_valid safe_line_count safe_grep_count
export LINDEMANN_SCRIPT DUMP2XYZ_SCRIPT SIMPLIFY_SCRIPT CONVERGENCE_SCRIPT UNWRAP_ASE_SCRIPT PER_ELEMENT_SCRIPT
export NEED_PREPROCESS DUMP2XYZ_ARGS SIMPLIFY_ARGS
export USE_UNWRAP UNWRAP_METHOD UNWRAP_ELEMENTS
export ATOM_SELECTION DO_COMPARISON SKIP_EXISTING CONVERGENCE_THRESHOLD
export ENABLE_PER_ELEMENT
export MASTER_CSV COMPARISON_CSV CONVERGENCE_CSV SUMMARY_LOG ERROR_LOG PROGRESS_FILE
export MASTER_LOCK COMPARISON_LOCK CONVERGENCE_LOCK PROGRESS_LOCK
export OUTPUT_MODE total

# ==================== 初始化 ====================

{
    echo "========================================"
    echo "纳米团簇Lindemann指数批量计算 v8"
    echo "支持分元素分析 (Pt, Sn, PtSn, PtSnO)"
    echo "========================================"
    echo "运行ID: $RUN_ID"
    echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "总目录数: $total"
    echo "并行数: $MAX_JOBS"
    echo ""
    echo "配置:"
    echo "  预处理: $NEED_PREPROCESS"
    echo "  Unwrap: $USE_UNWRAP ($UNWRAP_METHOD)"
    echo "  对比计算: $DO_COMPARISON"
    echo "  分元素计算: $ENABLE_PER_ELEMENT"
    echo "  跳过已有: $SKIP_EXISTING"
    echo ""
} > "$SUMMARY_LOG"

cat << EOF
========================================
纳米团簇Lindemann指数批量计算 v8
分元素分析: Pt, Sn, PtSn, PtSnO
========================================
运行ID: $RUN_ID
总目录数: $total
并行数: $MAX_JOBS

开始处理...
EOF

# ==================== 进度监控 ====================

if [ "$SHOW_PROGRESS" = true ]; then
    (
        while true; do
            sleep $UPDATE_INTERVAL

            current=$(cat "$PROGRESS_FILE" 2>/dev/null | head -1 | tr -d '[:space:]' || echo "0")
            if [[ ! "$current" =~ ^[0-9]+$ ]] || [ -z "$current" ]; then
                current=0
            fi

            if [ "$current" -gt 0 ] && [ "$total" -gt 0 ]; then
                percent=$(awk "BEGIN {printf \"%.1f\", $current/$total*100}" 2>/dev/null || echo "0.0")
            else
                percent="0.0"
            fi

            success=0
            skip=0
            if [ -f "$MASTER_CSV" ]; then
                total_lines=$(safe_line_count "$MASTER_CSV")
                if [ "$total_lines" -gt 1 ]; then
                    all_count=$((total_lines - 1))
                    skip=$(safe_grep_count ",skip," "$MASTER_CSV")
                    success=$((all_count - skip))

                    if [ "$success" -lt 0 ]; then success=0; fi
                    if [ "$skip" -lt 0 ]; then skip=0; fi
                fi
            fi

            fail=$((current - success - skip))
            if [ "$fail" -lt 0 ]; then fail=0; fi

            echo -ne "\r进度: $current/$total ($percent%) | ✅ $success | ⏭️  $skip | ❌ $fail    "

            if [ "$current" -ge "$total" ]; then
                break
            fi
        done
    ) &
    PROGRESS_PID=$!
fi

# ==================== 主循环 ====================

if [ "$USE_PARALLEL" = true ]; then
    if command -v parallel &> /dev/null; then
        printf "%s\n" "${directories[@]}" | \
            parallel -j "$MAX_JOBS" --line-buffer \
            process_single_directory {} {#} > /dev/null
    else
        for i in "${!directories[@]}"; do
            echo "${directories[$i]} $((i+1))"
        done | \
            xargs -P "$MAX_JOBS" -n 2 bash -c 'process_single_directory "$0" "$1"' > /dev/null
    fi
else
    for i in "${!directories[@]}"; do
        process_single_directory "${directories[$i]}" "$((i+1))"
    done
fi

if [ "$SHOW_PROGRESS" = true ] && [ -n "${PROGRESS_PID:-}" ]; then
    kill $PROGRESS_PID 2>/dev/null || true
    wait $PROGRESS_PID 2>/dev/null || true
fi

echo ""

# ==================== 统计 ====================

success_count=0
skip_count=0

if [ -f "$MASTER_CSV" ]; then
    total_lines=$(safe_line_count "$MASTER_CSV")

    if [ "$total_lines" -gt 1 ]; then
        all_count=$((total_lines - 1))
        skip_count=$(safe_grep_count ",skip," "$MASTER_CSV")
        success_count=$((all_count - skip_count))

        if [ "$success_count" -lt 0 ]; then success_count=0; fi
        if [ "$skip_count" -lt 0 ]; then skip_count=0; fi
    fi
fi

fail_count=$((total - success_count - skip_count))
if [ "$fail_count" -lt 0 ]; then fail_count=0; fi

if [ "$total" -gt 0 ]; then
    success_rate=$(awk "BEGIN {printf \"%.1f\", $success_count/$total*100}")
else
    success_rate="0.0"
fi

# ==================== 总结 ====================

cat << EOF

========================================
🎉 处理完成
========================================
完成时间: $(date '+%Y-%m-%d %H:%M:%S')

统计:
  总任务: $total
  ✅ 成功: $success_count
  ⏭️  跳过: $skip_count
  ❌ 失败: $fail_count
  成功率: ${success_rate}%

EOF

if [ "$DO_COMPARISON" = true ] && [ -f "$COMPARISON_CSV" ]; then
    comp_count=$(safe_line_count "$COMPARISON_CSV")

    if [ "$comp_count" -gt 1 ]; then
        echo "========================================="
        echo "分元素结果预览（前10条）"
        echo "========================================="

        tail -n +2 "$COMPARISON_CSV" 2>/dev/null | head -10 | \
            awk -F',' '{printf "%-40s T%-5sK | Cluster=%.6f | Pt=%.6f Sn=%.6f PtSn=%.6f PtSnO=%.6f\n", 
                $2, $3, $4, $8, $9, $10, $11}'

        echo ""

        # 计算各元素平均值
        avg_cluster=$(tail -n +2 "$COMPARISON_CSV" 2>/dev/null | \
            awk -F',' '{sum+=$4; n++} END {if(n>0) printf "%.6f", sum/n; else print "N/A"}')
        avg_pt=$(tail -n +2 "$COMPARISON_CSV" 2>/dev/null | \
            awk -F',' '{sum+=$8; n++} END {if(n>0) printf "%.6f", sum/n; else print "N/A"}')
        avg_sn=$(tail -n +2 "$COMPARISON_CSV" 2>/dev/null | \
            awk -F',' '{sum+=$9; n++} END {if(n>0) printf "%.6f", sum/n; else print "N/A"}')
        avg_ptsn=$(tail -n +2 "$COMPARISON_CSV" 2>/dev/null | \
            awk -F',' '{sum+=$10; n++} END {if(n>0) printf "%.6f", sum/n; else print "N/A"}')
        avg_ptsno=$(tail -n +2 "$COMPARISON_CSV" 2>/dev/null | \
            awk -F',' '{sum+=$11; n++} END {if(n>0) printf "%.6f", sum/n; else print "N/A"}')

        echo "平均Lindemann指数:"
        echo "  Cluster(原方法): $avg_cluster"
        echo "  Pt:             $avg_pt"
        echo "  Sn:             $avg_sn"
        echo "  PtSn:           $avg_ptsn"
        echo "  PtSnO:          $avg_ptsno"
        echo ""
    fi
fi

cat << EOF
输出文件:
  主表: $MASTER_CSV
  对比(含分元素): $COMPARISON_CSV
  收敛: $CONVERGENCE_CSV
  日志: $SUMMARY_LOG
  错误: $ERROR_LOG

查看命令:
  cat $COMPARISON_CSV
========================================
EOF

rmdir "$MASTER_LOCK" "$COMPARISON_LOCK" "$CONVERGENCE_LOCK" "$PROGRESS_LOCK" 2>/dev/null || true

exit 0
