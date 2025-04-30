source venv/bin/activate
cd Spare
pip install -e .
cd ..
./run_scripts/run_unqomp.sh 1
./run_scripts/run_hand.sh 1
./run_scripts/run_spire.sh 1
./run_scripts/run_length_simplified.sh 1
