BENCHMARKS=(
 twooffive
 #nam_nam_vbe_adder_3_check_check
 #nam_hlf_5 
 #nam_vbe_adder_3
 #nam_vbe_adder_3_check
 #nam_csla_mux_3
 ##nam_csum_mux_9_check
 #nam_csum_mux_9
 #nam_csla_mux_3_check
 #nam_mod5_4
 #nam_mod5_4_check
 #nam_mod_mult_55
 #nam_mod_mult_55_check_check
 # nam_ham15-low
 #nam_muliply_5
 # nam_hlf_5
 # nam_QFT8
 #nam_Adder16
 
 #nam_mod_mult_55
 # nam_mod_mult_55_check
 # nam_mod5_4_check
 #nam_gf2^4_mult
 #nam_gf2^4_mult_check
 #toffoli_3
 #toffoli_4
 #toffoli_5
 #toffoli_6
 #incrementer_4
 #incrementer_5
 #takahashi_4 
 #takahashi_6
 #takahashi_8 
 #cuccaro_4
 #cuccaro_6
 #random_circuit_6
 #random_circuit_11
 #random_circuit_16
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
	mkdir -p $(pwd)/"OUTPUT_JSON"/"TOOL_RESULTS_ANCILLA"
	dotpath=$DOT_PATH/${bench}/
	decomplogpath=$LOG_PATH/${bench}/"decomplog3.txt"
	rewritelogpath=$LOG_PATH/${bench}/"rewriteloge2.txt"
	json_path=$(pwd)/"OUTPUT_JSON"/"TOOL_RESULTS_ANCILLA"
	echo "Testing $bench..."
	echo "CHECK START"	

	#rm -r $dotpath/*.json
	#rm -r $dotpath/*.txt
	#rm -r $dotpath/*.dot

	#pytest $decomp/decomposition/decomp_tool_runner.py --dotpath $dotpath \
	#			--circuit-name $bench --benchmark-json=$json_path/output_$bench.json -s | tee $decomplogpath 	
	# status_decomp=0
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
			python $rewriter/compile_rewrite_graphs.py --test-name=$dotpath${bench}_graph --compile-mode="optimized"\
			 	--runs=$1 --json-name=$json_path/output_$bench.json >  $rewritelogpath
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
