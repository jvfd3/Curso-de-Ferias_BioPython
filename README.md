# Curso de Férias: Biopython

## Material do curso

- **Nosso conteúdo**
  - Apresentação
    - [Canva](https://www.canva.com/design/DAG_osc-YGg/sWj0-_oCf7b9Y76z4LHICQ/view)
    - [PDF](<../Documentos/Minicurso - Biopython.pdf>)
  - [GitHub](https://github.com/jvfd3/Curso-de-Ferias_BioPython)
- **Conteúdo base:**
  - [Documentação do biopython](https://biopython.org)
  - [Tutorial e Cookbook](https://biopython.org/docs/latest/Tutorial/index.html)
- **Conteúdo adicional no GitHub:**
  - [peterjc/biopython_workshop](https://github.com/peterjc/biopython_workshop/blob/master/reading_sequence_files/README.rst)
  - [tiagoantao/biopython-notebook](https://github.com/tiagoantao/biopython-notebook)
  - [peterjc/biopython_workshop](https://github.com/peterjc/biopython_workshop)

## Módulos utilizados

- Bio.Seq
- Bio.SeqUtils
- Bio.SeqIO
- Bio.SeqRecord
- Bio.SeqFeature
- Bio.Align
- Bio.PDB

## Sumário

- **1.** Preparo do material
  - **1.1.** Instalando Biopython
  - **1.2.** Obtendo FASTA do NCBI
  - **1.3.** Obtendo arquivo do ExPASy
  - **1.4.** Obtendo arquivos diretamente da internet
- **2.** Sequências
  - **2.1.** Objeto Seq
    - **2.1.1.** Métodos do Seq
    - **2.3.2.** Busca de Subsequências
    - **2.4.1.** Lendo arquivos
      - **2.4.1.1.** Lendo arquivo FASTA
      - **2.4.1.2.** Lendo arquivo GenBank
      - **2.4.1.3.** Lendo arquivos compactados (gzip)
      - **2.4.1.4.** Carregando o arquivo para a memória
      - **2.4.1.5.** Percorrendo registros
    - **2.4.2.** Escrevendo/Convertendo arquivos
    - **2.5.1.** Fatiando um SeqRecord
    - **2.5.2.** Escrevendo arquivos
  - **2.6.** Seq Feature
    - **2.6.1.** Exemplo real do Seq Feature em um arquivo GenBank
- **3.** Alinhamento de Sequências: Bio Align
  - **3.1.** Criando e inspecionando um Alignment
  - **3.2.** Alinhamento par-a-par com PairwiseAligner
  - **3.3.** Alinhamento de proteínas com matriz de substituição
  - **3.4.** Alinhamento múltiplo
  - **3.5.** Leitura de alinhamentos de arquivo
- **4.** Estrutural: Bio.PDB
  - **4.1.** Percorrendo átomos, resíduos e cadeias
  - **4.2.** Lendo e escrevendo arquivos estruturais
    - **4.2.1.** Escrevendo arquivos em CIF
    - **4.2.2.** Escrevendo arquivos em PDB
  - **4.3.** Análise estrutural (SMCRA)
    - **4.3.1.** Aprofundando nos átomos
      - **4.3.1.1.** Coordenadas internas: distâncias, ângulos, ângulos de torção
  - **4.4.** Modificando coordenadas internas de proteínas
    - **4.4.1.** Criando mapa de distâncias
  - **4.5.** Recriando mapa de distâncias
  - **4.6.** Impressão 3D de estruturas de proteínas
