FOLDER="benches_unqomp/"

BENCHMARKS=(
 adder_6_anc
 adder_12_anc
 integercomparator_6_anc
 integercomparator_12_anc
 wa_4_anc
 wa_12_anc
 mult_3
 mult_6
 mult_12
 wassaveqb_4
 wassaveqb_8
 wassaveqb_12
 mcry_12
 dj_10
 grover_8
 plr_3_anc
 plr_12_anc
)

# Vars
export PYTHONPATH=$(pwd)
export DOT_PATH="$(pwd)/dot_files/"
export LOG_PATH="$(pwd)/log_files/"

errors=()
RED='\033[0;31m'
NC='\033[0m' # No Color
export RIPL_HOME=$basedir
rewriter=$(pwd)/Spare/rewriter


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
	rewritelogpath=$LOG_PATH/$FOLDER/${bench}/"rewritelog7.txt"
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
