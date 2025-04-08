BENCHMARKS=(
 #adder2_anc
 adder6_anc
 #adder12_anc
 #intergercomparator_6_anc
 #intergercomparator_12_anc
 #wa_4_anc
 #wa_12_anc
 #mult_3
 #mult_6
 #mult_12
 #t_wassaveqb_4
 #t_wassaveqb_8
 #t_wassaveqb_12
 #mcry_12
 #temp_adder6
 #temp_adder12
 #temp_intergercomparator_6_20
 #temp_intergercomparator_12_40
 #temp_wa_4
 #temp_wa_12
 #temp_mult_3
 #temp_mult_6
 #temp_mult_12
 #mult_3
 #mult_6
 #mult_12
 #temp_wassaveqb_4
 #temp_wassaveqb_8
 #temp_wassaveqb_12
 # mcx_12 
 #dj_10
 #grover_8
 
 #t_pld_3_anc
 #t_pld_12_anc
 #new_pld_12_anc
 #new_pld_3_anc
 # temp_wassaveqb_4
 
 #t_wassaveqb_2
 #t_wassaveqb_8
 #t_wassaveqb_12
 #wassaveqb_4
 #wassaveqb_8
 #wassaveqb_12
 #length_simplified_orig1
 #length_simplified1
 #length_simplified_only_cn1
 #length_simplified_only_cf1
 #length_simplified5
 #length_simplified_only_cn5
 #length_simplified_only_cf5
 #length_simplified_orig5
 #length_simplified9
 #length_simplified_orig9
 #length_simplified_only_cf9
 #length_simplified_only_cn9
 #simpl_find_pos1
 #simpl_find_pos_orig1
 #simpl_is_prefix1
 #simpl_is_prefix_orig1
 #simpl_pop_front1 #length1
 #simpl_length1
 #simpl_push_back1
 #simpl_push_back_orig1
 #simpl_sum1
 #simpl_compare1
 #simpl_remove1
 #simpl_num_matching1
 ###########
 #twooffive
 #rd53
 #sixsim
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
	mkdir -p $(pwd)/"OUTPUT_JSON"/"TOOL_RESULTS_Unqomp2"
	dotpath=$DOT_PATH/${bench}/
	decomplogpath=$LOG_PATH/${bench}/"decomplog3.txt"
	rewritelogpath=$LOG_PATH/${bench}/"rewritelog_opt_check.txt"
	json_path=$(pwd)/"OUTPUT_JSON"/"TOOL_RESULTS_Unqomp_time_only"
	echo "Testing $bench..."
 	
	#rm $dotpath/*.txt
	#rm $dotpath/*.dot
	#pytest $decomp/decomposition/decomp_tool_runner.py --dotpath $dotpath \
	#	--circuit-name $bench --return-base --benchmark-json=$json_path/output_$bench.json -s > $decomplogpath 	
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
