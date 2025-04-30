BENCHMARKS=(
# new_adder6
# new_adder12
# # adder12_anc
adder6_anc
# new_intergercomparator_12 #_anc
# new_intergercomparator_6 # _anc
# new_wa_12 # _anc
# new_wa_4 #_anc
#  wassaveqb_12
wassaveqb_4
#  wassaveqb_8
#  new_wassaveqb_12
#  new_wassaveqb_8
#  new_wassaveqb_4
# intergercomparator_12_anc
# wa_12
# wa_4
# new_mult_6
# new_mult_3
# new_mult_12
mult_6
# mult_3
# mult_12
#comp_length1
#comp_sum1
#comp_find_pos1
#comp_pop_front1
#comp_remove1
#simpl_length1
#simpl_insert1
#simpl_contains1
# cp_length1
#new_pld_3_anc
#new_pld_12_anc
#_sum9
#comp_num_matching1
#simpl_find_pos9
#simpl_is_prefix9
#simpl_push_back9
#simpl_compare9
#comp_compare1
# length_simplified_orig1
# length_simplified1
# length_simplified_only_cn1
#length_simplified_only_cf1
#  length_simplified5
#  length_simplified_only_cn5
#  length_simplified_only_cf5
#  length_simplified_orig5
#  length_simplified_orig9
#  length_simplified_only_cf9
#  length_simplified_only_cn9
#  length_simplified9
#  length1
#  length5
#  length9
#  length_orig1
#  length_orig5
#  length_orig9
#  compare5
#  compare9
#  compare_orig1
#  compare_orig5
#  compare_orig9
# find_pos5
# find_pos9
# find_pos_orig1
#  find_pos_orig5
#  find_pos_orig9
#  is_prefix5
# is_prefix9
#comp_is_prefix1 # _orig1
#  is_prefix_orig5
#  is_prefix_orig9
#comp_push_back1
# push_back5
# push_back9
# push_back_orig1
#  push_back_orig5
#  push_back_orig9
###############
#errorcorr
#adder6
#twooffive
#sixsim
#rd53
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
mkdir -p $DOT_PATH/"manual_circs"
mkdir -p $LOG_PATH

for b in ${!BENCHMARKS[@]}; do
	bench=${BENCHMARKS[$b]}
	mkdir -p $DOT_PATH/${bench}
	mkdir -p $LOG_PATH/${bench}
	dotpath=$DOT_PATH/"manual_circs"/
	decomplogpath=$LOG_PATH/${bench}/"manual_decomplog_e2.txt"
	rewritelogpath=$LOG_PATH/${bench}/"manual_simulation_log.txt"
	json_path=$(pwd)/"OUTPUT_JSON"/"MANUAL_RESULTS_IMPROVED"
	echo "Testing $bench..."
	echo $rewritelogpath
	if [$2 -eq 2]; then
		python $rewriter/compile_rewrite_graphs.py --test-name=$dotpath${bench}_graph \
			--compile-mode="bench_qubit" --runs=$1 --dont-final-replace --json-name=$json_path/output_$bench.json > $rewritelogpath
	else
		python $rewriter/compile_rewrite_graphs.py --test-name=$dotpath${bench}_graph \
			--compile-mode="bench_qubit" --runs=$1 --json-name=$json_path/output_$bench.json > $rewritelogpath
	fi
	status=$?
	if [ $status -gt 0 ]
	then 
	  errors+=("${line}, ${bench}")
	fi
done

echo -e "${RED}Failed tests:"
for i in ${!errors[@]}; do
    error=${errors[$i]} 
    echo -e "${RED}$error,"
done
echo -e "${NC}"
