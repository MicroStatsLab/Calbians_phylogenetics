# Normalize VCF to biallelic
bcftools norm -m -any 250215_All_738_high_quality_vcf_SNPs.vcf.gz -Oz -o vcf_biallelic.vcf.gz
bcftools index vcf_biallelic.vcf.gz

# Convert PL to GL, handling missing or zero-sum values
bcftools query -f '%CHROM:%POS\t%REF\t%ALT[\t%PL]\n' vcf_biallelic.vcf.gz \
| perl -F'\t' -ane '
    @out = ($F[0], $F[1], $F[2]);
    for ($i = 3; $i < @F; $i++) {
        @pl = split(",", $F[$i]);
        if (@pl != 3) { @pl = (".", ".", "."); }
        @gl = map { ($_ eq "." ? 0 : 10 ** (-$_/10)) } @pl;
        $sum = 0; $sum += $_ for @gl;
        if ($sum == 0) { @gl = (1/3, 1/3, 1/3); }
        push @out, @gl;
    }
    print join("\t", @out), "\n";
' > output.beagle

# Generate header
bcftools query -l vcf_biallelic.vcf.gz > samples.txt
(
    printf "marker\tallele1\tallele2\n"
    cat samples.txt | while read s; do
        printf "%s\t%s\t%s\t" "$s" "$s" "$s"
    done
) | sed 's/\t$//' > header.txt

# Concatenate header and GLs
cat header.txt output.beagle > final.beagle

# Verify
awk -F'\t' '{print NF; exit}' final.beagle

