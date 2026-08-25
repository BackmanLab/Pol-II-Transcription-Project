#!/bin/bash
##################################################################
# Author: Lucas Carter
# PI: Vadim Backman
# Description:
# Quest-specific adaptation of mega.sh. Run once from the Quest
# login node: bash patch_mega.sh
#
# Applies:
#   - default -A account flag
#   - increases HIC job memory to 400G
#   - makes HIC job thread count follow the -T flag
#
# Figure output: none -- infrastructure/patch script.
##################################################################

MEGA=/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer/scripts/mega.sh

python3 - << 'PYEOF'
path = "/projects/b1042/BackmanLab/Lucas/060825_hic-samp/juicer/scripts/mega.sh"
with open(path) as f:
    c = f.read()

# 1. Add default account variable after groupname line
c = c.replace(
    'groupname="a$(date +%s)"',
    'groupname="a$(date +%s)"\naccount="b1042"',
    1
)

# 2. Add -A to getopts string
c = c.replace(
    '"d:g:hfs:S:l:L:q:Q:D:y:T:"',
    '"d:g:hfs:S:l:L:q:Q:D:y:T:A:"',
    1
)

# 3. Add account case handler
c = c.replace(
    'T) threadsHic=$OPTARG ;;',
    'T) threadsHic=$OPTARG ;;\n        A) account=$OPTARG ;;',
    1
)

# 4. Replace all hardcoded --account=b1042 with variable
c = c.replace('--account=b1042', '--account=${account}')

# 5. HIC1 and HIC30 jobs: increase memory 150G -> 400G
c = c.replace(
    'J "${groupname}_hic1"\n#SBATCH -d "${dependinter1}"\n#SBATCH --mem=150G',
    'J "${groupname}_hic1"\n#SBATCH -d "${dependinter1}"\n#SBATCH --mem=400G'
)
c = c.replace(
    'J "${groupname}_hic30"\n#SBATCH -d "${dependinter30}"\n#SBATCH --mem=400G',
    'J "${groupname}_hic30"\n#SBATCH -d "${dependinter30}"\n#SBATCH --mem=400G'
)

# 6. HIC1 and HIC30 jobs: make -c follow threadsHic
c = c.replace(
    '#SBATCH -c 8\n#SBATCH --ntasks=1\n#SBATCH -o $logdir/hic1',
    '#SBATCH -c ${threadsHic:-8}\n#SBATCH --ntasks=1\n#SBATCH -o $logdir/hic1'
)
c = c.replace(
    '#SBATCH -c 8\n#SBATCH --ntasks=1\n#SBATCH -o $logdir/hic30',
    '#SBATCH -c ${threadsHic:-8}\n#SBATCH --ntasks=1\n#SBATCH -o $logdir/hic30'
)

# 7. Update Java heap to 400G in hic jobs
c = c.replace(
    'export IBM_JAVA_OPTIONS="-Xmx150000m -Xgcthreads1"\nexport _JAVA_OPTIONS="-Xms150000m -Xmx150000m"',
    'export IBM_JAVA_OPTIONS="-Xmx400000m -Xgcthreads1"\nexport _JAVA_OPTIONS="-Xms400000m -Xmx400000m"'
)

with open(path, 'w') as f:
    f.write(c)
print("Done")
PYEOF

# Verify all patches applied
echo "=== account default variable ==="
grep -n "^account=" $MEGA

echo "=== getopts string ==="
grep -n "getopts" $MEGA

echo "=== -A case handler ==="
grep -n "account=\$OPTARG" $MEGA

echo "=== sbatch account (first 3) ==="
grep -n "\-\-account=" $MEGA | head -3

echo "=== HIC job memory and threads ==="
grep -n "hic1\|hic30\|SBATCH -c\|mem=400\|Xmx400" $MEGA | head -15

echo ""
bash -n $MEGA && echo "Syntax OK"
