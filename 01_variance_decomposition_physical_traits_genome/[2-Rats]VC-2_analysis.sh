
### STEP2. Run GCTA ###

gcta=/SAN/mottlab/heterosis/bin/gcta/gcta64
plink=/home/ucbtlcu/bin/plink2
wd=/SAN/mottlab/heterosis/3_HSrats/5VC/AllInds
datdir=/SAN/mottlab/heterosis/3_HSrats/5VC/data

cd $wd

traits=(`cat $datdir/HSrats_res.traitnames`)
let N=${#traits[@]}-1 
idx=(`seq 0 $N`)

$plink --noweb --tfile $datdir/HSrats_res --maf 0.05 --geno 0.1 --make-bed --out HSrats_Final

$gcta --bfile HSrats_Final --thread-num 10  --autosome --make-grm-gz --out HSrats_Final_add
$gcta --bfile HSrats_Final --thread-num 10  --autosome --make-grm-d-gz --out HSrats_Final_dom

$gcta --bfile HSrats_Final --thread-num 10  --autosome --make-grm --out HSrats_Final_add
$gcta --bfile HSrats_Final --thread-num 10  --autosome --make-grm-d --out HSrats_Final_dom

mkdir add_dom_VC

for i in ${idx[@]};do
    let j=$i+1
    $gcta --reml --mgrm $datdir/add_domi.txt --pheno $datdir/HSrats_res.pheno --mpheno $j --out ./add_dom_VC/${traits[$i]}
done


