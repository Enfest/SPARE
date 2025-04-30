FOLDER="benches_unqomp/"

BENCHMARKS=(
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
export PYTHONPATH=$(pwd)
export DOT_PATH="$(pwd)/dot_files/"
export LOG_PATH="$(pwd)/log_files/"

errors=()
RED='\033[0;31m'
NC='\033[0m' # No Color

basedir=$(pwd)
rewriter=$(pwd)/Spare/rewriter
export RIPL_HOME=$basedir


mkdir -p "$DOT_PATH"
mkdir -p $LOG_PATH

for b in ${!BENCHMARKS[@]}; do
	bench=${BENCHMARKS[$b]}
	mkdir -p $DOT_PATH/$FOLDER/${bench}
	mkdir -p $LOG_PATH/$FOLDER/${bench}
	mkdir -p $(pwd)/"OUTPUT_JSON"/
	mkdir -p $(pwd)/"OUTPUT_JSON"/"TOOL_RESULTS"
	mkdir -p $(pwd)/"temp_files"/
	dotpath=$DOT_PATH/$FOLDER/${bench}/
	rewritelogpath=$LOG_PATH/$FOLDER/${bench}/"rewritelog.txt"
	json_path=$(pwd)/"OUTPUT_JSON"/"TOOL_RESULTS"
	echo "Testing $bench..."
	echo $rewriter/compile_rewrite_graphs.py
	python $rewriter/compile_rewrite_graphs.py --test-name=$dotpath${bench}_graph --compile-mode="optimized_qubit"\
		--runs=$1 --json-name=$json_path/output_$bench.json > $rewritelogpath
	rm -r $(pwd)/"temp_files"/
done

# python $basedir/src/bench_csv_aggregator.py $ $basedir/$benchout/suitesparse_check_sweep_reuse_${2}_${bench}.csv
echo -e "${RED}Failed tests:"
for i in ${!errors[@]}; do
    error=${errors[$i]} 
    echo -e "${RED}$error,"
done
echo -e "${NC}"
