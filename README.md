# Anotações

## Plano de Aula: Biopython

### Introdução a Biopython

27/09/2025 - 12:47:10

[milennapirovani@hotmail.com](mailto:milennapirovani@hotmail.com)

Milenna Machado Pirovani

PPG Bioinformática

Teórico e Prático: 8 horas/aula (duas manhãs, de 09 às 12:20hrs)

Sim

Sim, o minicurso será totalmente prático e será fundamental que os participantes acompanhem as atividades em sala

Capacidade máxima da sala disponibilizada para o minicurso.

Sim

Python básico (listas, dicionários, funções), Jupyter Notebook/Colab ou terminal; laptop com Python 3.8+; conda/pip instalados.

Sim, de 1 a 3 monitores

#### Objetivos do Minicurso

Capacitar os participantes a usar Biopython para tarefas práticas (aplicado a genômica, transcriptômica e proteômica, apresentando as funções mais básicas da biblioteca).

Construir pipelines simples reprodutíveis (download → processamento → análise → relatório).

Desenvolver autonomia para integrar Biopython com ferramentas externas (MSA, BLAST, aligners) e gerar outputs interpretáveis (FASTA/GenBank/MSA/árvores/PDB parsing).

#### Conteúdo programático

##### Módulo 1: Genoma e Transcriptoma (Dia 1)

- **Conceitos**: Seq, SeqRecord, SeqIO, formatos (FASTA vs GenBank);
- **Parse de GenBank**: _SeqRecord.features_, extração de CDS/genes, metadados;
- **Análises básicas**: GC%, contagem de k-mers, estatísticas de comprimento;
- **Exercício prático**: 1. Baixar um GenBank; 2. Extrair CDS; 3. Traduzir; 4. Salvar FASTA;
- **FASTQ**: leitura com SeqIO, scores de qualidade, conversão FASTQ→FASTA;
- Filtragem de reads por qualidade/comprimento; resumo estatístico de reads;
- **Alinhamento Múltiplo (MSA)**: AlignIO, usar MAFFT/Clustal (wrapper ou subprocess), extrair consenso;
- Construção básica de árvore filogenética a partir de MSA (uso de Bio.Phylo para visualização);
- **Exercício prático**: pipeline FASTQ → filtragem → MSA → análise de consenso.

##### Módulo 2: Proteoma (Dia 2, tarde)

- **Tradução de CDS:** proteínas; criação e manipulação de SeqRecord proteicos;
- _Bio.SeqUtils.ProtParam.ProteinAnalysis_: massa molecular, pI estimado, composição aminoacídica, índices;
- _Bio.SearchIO_: parsing de resultados BLAST / HMMER / Search outputs.
- Bio.PDB: leitura de PDB, acesso a cadeias/resíduos, cálculos de distância/RMSD (conceito e exemplo simples);
- **Exercício prático**: análise proteica (massa/pI/composição) + buscar e analisar PDB relacionado.

#### Metodologia e formato das aulas

- Cada bloco:
  - Teoria curta (10-15 min);
  - Apresentação ao vivo em Jupyter (15-25 min);
  - Exercício prático (20-40 min).
- Notebooks com células "execute aqui" e gabaritos; datasets de exemplo (pequenos) distribuídos;
- Trabalhos em pares ou individuais durante o minicurso;
- Material pré-curso: checklist de instalação (conda/pip), links para instalação de MAFFT/Clustal/BLAST (opcional).
