BENCHMARKS=(
 # toffoli_3
 # incrementer_5 
 # toffoli_5
 # toffoli_6
 #incrementer_4
 #incrementer_5
 # takahashi_4
 takahashi_6
 # benchmarks still being seen
 # takahashi_8
 # cuccaro_4
 # cuccaro_6
 # random_circuit_6
 # random_circuit_11
 # random_circuit_16
)

# Vars
export PYTHON_PATH=$(pwd)
export DOT_PATH="$(pwd)/dot_files/"
export LOG_PATH="$(pwd)/log_files/"

errors=()
RED='\033[0;31m'
NC='\033[0m' # No Color

# format_outdir=${SUITESPARSE_FORMATTED_PATH} 
basedir=$(pwd)
decomp=$(pwd)/decomp
rewriter=$(pwd)/rewriter
export RIPL_HOME=$basedir


mkdir -p "$DOT_PATH"
mkdir -p $LOG_PATH

for b in ${!BENCHMARKS[@]}; do
	bench=${BENCHMARKS[$b]}
	mkdir -p $DOT_PATH/${bench}
	mkdir -p $LOG_PATH/${bench}
	mkdir -p $(pwd)/"OUTPUT_JSON"/
	mkdir -p $(pwd)/"OUTPUT_JSON"/"TOOL_RESULTS"
	dotpath=$DOT_PATH/${bench}/
	decomplogpath=$LOG_PATH/${bench}/"decomplog3.txt"
	rewritelogpath=$LOG_PATH/${bench}/"rewritelog2.txt"
	json_path=$(pwd)/"OUTPUT_JSON"/"TOOL_RESULTS"
	echo "Testing $bench..."

 	# pytest $decomp/test/test_decomp_tool.py --dotpath $dotpath --circuit-name $bench --benchmark-json=$json_path/output_$bench.json -s > $decomplogpath 
	if [ $2 -eq 2 ]; then
		if [$3 -eq 2]; then
			python $rewriter/compile_rewrite_graphs --test-name=${bench}_graph --compile-mode="bench" --runs=$1 --dont-final-replace > $logpath
		else
			python $rewriter/compile_rewrite_graphs --test-name=${bench}_graph --compile-mode="bench" --runs=$1 > $logpath
		fi
	else
		printf $rewriter/compile_rewrite_graphs.py "\n"
		python $rewriter/compile_rewrite_graphs.py --test-name=$dotpath${bench}_graph --compile-mode="optimized" --runs=$1 --json-name=$json_path/output_$bench.json > $rewritelogpath
	fi
	#python $basedir/scripts/converter.py --json_name $path/$line.json
	status=$?
	if [ $status -gt 0 ]
	then 
	  errors+=("${line}, ${bench}")
	fi
done

# python $basedir/src/bench_csv_aggregator.py $ $basedir/$benchout/suitesparse_check_sweep_reuse_${2}_${bench}.csv
echo -e "${RED}Failed tests:"
for i in ${!errors[@]}; do
    error=${errors[$i]} 
    echo -e "${RED}$error,"
done
echo -e "${NC}"
