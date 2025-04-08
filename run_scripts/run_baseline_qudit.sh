BENCHMARKS=(
#toffoli_3
#toffoli_4 
toffoli_5
toffoli_6
incrementer_4
incrementer_5
takahashi_4
takahashi_6
takahashi_8
cuccaro_4
cuccaro_6
random_circuit_6
random_circuit_11
random_circuit_16
nam_vbe_adder_3
nam_csla_mux_3
nam_csum_mux_9
nam_mod5_4
nam_mod_mult_55
nam_gf2^4_mult
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
	rewritelogpath=$LOG_PATH/${bench}/"manual_simulation_log_e2.txt"
	json_path=$(pwd)/"OUTPUT_JSON"/"MANUAL_RESULTS_IMPROVED"
	echo "Testing $bench..."
	# decomp/decomposition/qudit_decomposition_runner
	pytest $decomp/decomposition/qudit_decomposition_runner.py --dotpath $dotpath --circuit-name $bench  --benchmark-json=$json_path/output_manual$bench.json -s > $decomplogpath 
	if [ $1 -eq 2 ]; then
		if [$2 -eq 2]; then
			python $rewriter/compile_rewrite_graphs.py --test-name=$dotpath${bench}_graph \
				--compile-mode="bench" --runs=$1 --dont-final-replace --json-name=$json_path/output_$bench.json > $rewritelogpath
		else
			python $rewriter/compile_rewrite_graphs.py --test-name=$dotpath${bench}_graph \
				--compile-mode="bench" --runs=$1 --json-name=$json_path/output_$bench.json > $rewritelogpath
		fi
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
