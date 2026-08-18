InterProScan
Descrição

InterProScan
 é uma ferramenta de anotação funcional que escaneia sequências proteicas contra múltiplos bancos de domínios e famílias (Pfam, PANTHER, CDD, SMART, TIGRFAM, etc.).
Ele atribui funções prováveis, termos GO e vias metabólicas às proteínas preditas a partir do genoma.

No pipeline, o módulo consome proteínas geradas pelo GFFREAD e retorna arquivos em formatos padronizados (TSV e GFF3).

Requisitos

Container: interpro/interproscan:5.75-106.0

Banco de dados: baixado e descompactado previamente em
/temporario2/14191848/fork/BenchAnnot/data/eukaryotes/db/interproscan/data

Bind no container:
--bind /caminho/local/data:/opt/interproscan/data

Entrada

sample_id.faa (arquivo FASTA de proteínas preditas pelo GFFREAD)

Saída

sample_id.interpro.tsv → tabela com hits de domínios, GO e pathways.

sample_id.interpro.gff3 → anotação em formato GFF3 para integração com genoma.

Parâmetros principais

--ips_formats → formatos de saída (default: tsv,gff3)

--ips_test true|false

true: usa apenas as 10 primeiras sequências (teste rápido)

false: roda no proteoma completo

CPUs e memória definidos no nextflow.config:

cpus = 8

memory = '32 GB'

time = '12h' (ajuste conforme o cluster)

Execução no Nextflow
Teste rápido
nextflow run main.nf -entry gff_interpro --ips_test true

Execução completa
nextflow run main.nf --ips_test false

Comando interno

O módulo executa:

interproscan.sh \
  -i input.faa \
  -f tsv,gff3 \
  -cpu ${task.cpus} \
  -dp \
  -goterms \
  -pa \
  -b sample_id.interpro

Flags importantes

-cpu ${task.cpus} → usa o número de núcleos reservado pelo SLURM.

-dp → desabilita lookup online de matches pré-calculados (necessário em HPC sem internet).

-goterms → adiciona anotações GO.

-pa → adiciona vias metabólicas (pathways).

-b → define o prefixo base dos arquivos de saída.

Atenção: -b, -d e -o são exclusivos. Aqui usamos apenas -b.

Notas

Phobius, SignalP e TMHMM não estão habilitados (licença restrita).

Se quiser evitar avisos sobre eles nos logs, use:
-exclappl SignalP,TMHMM,Phobius

Crie o diretório de logs antes de rodar:

mkdir -p logs results/eukaryotes/interproscan