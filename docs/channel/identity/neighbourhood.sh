#!/bin/bash
# Genes around the channel's locus, from NCBI's public records (E-utilities), as TSV:
# start, end, strand, locus tag, GenBank product, RefSeq (WP_) accession and its name.
#
#   bash neighbourhood.sh > neighbourhood_CP013127.tsv
#
# Defaults: BD122_14570 on CP013127.1 (B. diazoefficiens USDA 122), +-20 kb.
ACC=${ACC:-CP013127.1}
FROM=${FROM:-3165000}
TO=${TO:-3207000}
E=https://eutils.ncbi.nlm.nih.gov/entrez/eutils

printf 'start\tend\tstrand\tlocus_tag\tgenbank_product\tprotein_id\trefseq\trefseq_name\n'
curl -s "$E/efetch.fcgi?db=nuccore&id=$ACC&rettype=fasta_cds_aa&retmode=text&seq_start=$FROM&seq_stop=$TO" |
grep '^>' | while read -r h; do
    lt=$(sed -n 's/.*\[locus_tag=\([^]]*\)\].*/\1/p' <<< "$h")
    prod=$(sed -n 's/.*\[protein=\([^]]*\)\].*/\1/p' <<< "$h")
    pid=$(sed -n 's/.*\[protein_id=\([^]]*\)\].*/\1/p' <<< "$h")
    loc=$(sed -n 's/.*\[location=\([^]]*\)\].*/\1/p' <<< "$h")
    strand=+; [[ $loc == complement* ]] && strand=-
    read -r a b <<< "$(tr -dc '0-9.' <<< "$loc" | sed 's/\.\./ /')"
    wp=""; name=""
    if [ -n "$pid" ]; then
        wp=$(curl -s "$E/efetch.fcgi?db=ipg&id=$pid&rettype=ipg&retmode=text" | awk -F'\t' '$7 ~ /^WP_/{print $7; exit}')
        [ -n "$wp" ] && name=$(curl -s "$E/efetch.fcgi?db=protein&id=$wp&rettype=gp&retmode=text" |
            awk '/^DEFINITION/{sub(/^DEFINITION +/,""); d=$0; next} /^ACCESSION/{print d; exit} d{d=d" "$0}' |
            tr -s ' ' | sed 's/^MULTISPECIES: //; s/ *$//; s/ \[[^]]*\]\.$//; s/\.$//')
        sleep 0.4                       # NCBI asks for <= 3 requests/s without an API key
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$a" "$b" "$strand" "$lt" "$prod" "$pid" "$wp" "$name"
done
