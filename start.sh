#!/bin/bash
echo ""
echo ""
echo -e "\033[0;32m====================\033[0m"
echo -e "\033[0;32mWelcome to FALAPhyl!\033[0m"
echo -e "\033[0;32m====================\033[0m"
echo ""
echo ""
echo "Copying the template files in your directly if they don't exist (environments.txt, input.yaml)'"
# Step 1: Copy files as needed (if applicable)
cp -p environments.txt data/
cp -p input.yaml data/

echo ""
echo ""
echo "Example command(s):"
echo -e "\033[0;36msnakemake alpha beta breakdown diff network subject_beta subject_alpha subject_diff --use-conda --cores all --keep-going --retries 5 --rerun-incomplete\033[0m"
echo "Please copy and paste the above command(s) or type your own command below."
echo ""
echo ""
# Open an interactive bash shell so the user can input commands
exec bash
