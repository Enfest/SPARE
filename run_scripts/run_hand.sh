BENCHMARKS=(
 #twooffive
 #sixsim
 errorcorr
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
	mkdir -p $(pwd)/"OUTPUT_JSON"/"TOOL_RESULTS_ANCILLA_HAND"
	dotpath=$DOT_PATH/${bench}/
	decomplogpath=$LOG_PATH/${bench}/"decomplog2.txt"
	rewritelogpath=$LOG_PATH/${bench}/"rewritelog2.txt"
	json_path=$(pwd)/"OUTPUT_JSON"/"TOOL_RESULTS_ANCILLA_HAND"
	echo "Testing $bench..."
 	
	#rm $dotpath/*.txt
	#rm $dotpath/*.dot
	pytest $decomp/decomposition/decomp_tool_runner.py --dotpath $dotpath \
		--circuit-name $bench --return-base --benchmark-json=$json_path/output_$bench.json -s > $decomplogpath 	
	status_decomp=0
	status_decomp=$?

	if [ $status_decomp -eq 0 ]
	    then
		echo "Decomp done"
	    if [ $2 -eq 2 ]
			then
			if [ $3 -eq 2 ]
			then
				python $rewriter/compile_rewrite_graphs --test-name=${bench}_graph \
					--compile-mode="bench" --runs=$1 --dont-final-replace > $logpath
			else
				python $rewriter/compile_rewrite_graphs --test-name=${bench}_graph \
					--compile-mode="bench" --runs=$1 > $logpath
			fi
	    else
			echo $rewriter/compile_rewrite_graphs.py
			python $rewriter/compile_rewrite_graphs.py --test-name=$dotpath${bench}_graph --compile-mode="optimized_qubit"\
			 	--runs=$1 --json-name=$json_path/output_$bench.json > $rewritelogpath
			status=$?
	   		if [ $status -gt 0 ]
	   			then
					errors+=("${line}, ${bench}")
	   			else
					echo "Rewrites done and written"
	   		fi
		fi
	else
		echo "Decomp failed"
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
