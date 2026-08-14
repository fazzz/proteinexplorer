# ProteinExplorer チュートリアル(日本語版)

`prot`の全コマンド群を実データで一通り触ってみます。使うのは**1A8O**
(HIV-1キャプシドC末端ドメイン、X線1.70Å、chain A、タンパク質残基70+水88
+セレノメチオニン)です。以下の出力はすべて実際にコマンドを実行して得た
ものをそのまま貼っています(手書きのサンプル出力ではありません)。

## データの取得

```bash
mkdir -p examples/1a8o
curl -sL -o examples/1a8o/1A8O.pdb \
  https://raw.githubusercontent.com/biopython/biopython/master/Tests/PDB/1A8O.pdb
```

1A8OはBiopython公式の`Bio.PDB`チュートリアルでも使われている、配布条件の
はっきりした実データです。

金属結合部位のデモなど、1A8Oには無い要素(金属イオンなど)を使う一部の
デモだけは`examples/illustrative/`に置いた小さな合成データを使います。
その都度「合成データです」と明記します。

## 1. import / info / status

```console
$ prot import examples/1a8o/1A8O.pdb --name 1a8o
Imported '1a8o' as p_bc8a53a4
  1 chain(s), 158 residue(s), 644 atom(s)
  hetero groups: MSE

$ prot info 1a8o
p_bc8a53a4  (1a8o)
  source: examples/1a8o/1A8O.pdb
  format: pdb   imported: 2026-08-13T23:48:41.549157+00:00
  method: x-ray diffraction
  resolution: 1.7 A
  models: 1   chains: 1
  residues: 158   atoms: 644
    protein: 70
    water: 88
  hetero groups: MSE
```

`import`のたびに`.proteinexplorer/`プロジェクトディレクトリが(無ければ)
作られ、ファイルはそのままコピーされます。`prot status`でこれまで
importした全構造体を一覧できます。

## 2. descriptor(記述子)

```console
$ prot descriptor 1a8o
p_bc8a53a4  (1a8o)
  molecular weight: 9016.3 Da (heavy-atom)
  atoms: 644   residues: 158   chains: 1
  ligands: 0   waters: 88
  SASA: 5451.8 A^2
  radius of gyration: 11.78 A
  contact density (CA-CA < 8A / residue): 4.17
  hydrophobic ratio: 0.43
  disulfide bonds: 1
  secondary structure (geometric): C=0.23, E=0.17, H=0.60
```

ジスルフィド結合が1本、純粋に幾何学的な判定(SG-SG距離)で見つかって
います。後段の`prot contact disulfide`でも同じ結合が出てきます。

## 3. 共通selection language

すべての解析コマンドは以下のようなselection式を受け付けます。

```
protein                          全タンパク質残基
chain A                          チェーン指定
resid 190:210                    残基番号の範囲
resname CYS                      残基名指定
atom CA                          原子名指定
chain A and backbone             論理演算の組み合わせ
within 5 (resname MSE)           ある選択範囲の近傍
```

## 4. geometry(幾何解析)

```console
$ prot geometry backbone-torsions 1a8o --chain A --resid 175
A/GLU175
  phi=-104.2  psi=4.4  omega=-179.2
  chi1=-55.2  chi2=-61.3  chi3=-37.5

$ prot geometry distance 1a8o "chain A and resid 155 and atom CA" "chain A and resid 220 and atom CA"
13.263 A
```

`prot geometry coords`(centroid/COM/bounding box/回転半径/plane fit/
principal axes/moment of inertia)、`prot geometry rmsd`、
`prot geometry distmatrix`も同じ構造データに対して使えます。

## 5. contact(接触解析)

```console
$ prot contact disulfide 1a8o
A/CYS198  --  A/CYS218   2.04 A
```

`prot contact hbond`/`saltbridge`/`hydrophobic`/`pipi`/`cationpi`/
`map`(残基間contact map)/`network`(全interaction typeを統合したedge
list)も同じ枠組みで使えます。

## 6. secondary(二次構造)

この実行環境にはDSSPが入っていないので、`--method auto`は依存なしの
phi/psi分類器に自動フォールバックします。

```console
$ prot secondary 1a8o --method geometric
p_bc8a53a4  (1a8o)   method=geometric
  A: CCEEEECCEEHHHHHHHHHHHHHHCCEEHHHHHHHHHCHHHHCEEHHHHHHHHCCCCCEEHHHHHHHCCC
  composition: C=0.23, E=0.17, H=0.60
```

外部ツールなしでもヘリックスに富んだ本来のfoldが再現できています。

## 7. plot(静的プロット)

```bash
prot plot ramachandran 1a8o rama.png
prot plot contact-map 1a8o contacts.png --mode heavy --cutoff 6
prot plot secondary 1a8o ss.png
```

![Ramachandranプロット](../examples/1a8o/outputs/ramachandran.png)

実データらしく、alpha領域とbeta領域にきれいに分布が集中しています。

![Contact map](../examples/1a8o/outputs/contact_map.png)

![Secondary structure](../examples/1a8o/outputs/secondary.png)

## 8. pocket(ポケット検出)

大きい構造では`--selection`で探索範囲を絞るのがおすすめです(グリッドが
爆発しないよう)。

```console
$ prot pocket detect 1a8o --selection "chain A and resid 190:210" --spacing 1.5
4 pocket(s) found:
  #1  volume=61 A^3  surface~=126 A^2  hydrophobicity=0.67  druggability~=0.69
      residues: A/ALA209, A/GLU187, A/GLY206, A/GLY208, A/LEU202, A/LYS203, A/MSE214, A/PRO207, A/VAL191
  #2  volume=37 A^3  surface~=99 A^2  hydrophobicity=0.50  druggability~=0.61
      residues: A/ASN195, A/ASP197, A/CYS198, A/CYS218, A/GLN155, A/GLU159, A/LYS158, A/PHE161, A/PRO157, A/PRO160
  ...
```

これは依存なしのLIGSITE風近似で、fpocket相当の本格実装ではありません。
druggabilityスコアも学習済みモデルではなく簡易ヒューリスティックです。

## 9. mutate(点変異)

```console
$ prot mutate 1a8o --chain A --resid 200 --to ALA
A/THR200 -> ALA
  method: cb_only
  atoms placed: N, CA, C, O, CB
  built-in fallback: backbone kept, C-beta idealized from N/CA/C. Full side
  chain beyond C-beta not built -- install Scwrl4 for a complete rotamer.
Saved as '1a8o_A200THRALA' (p_7d7a2c8f)
```

この環境にScwrl4が無いため、依存なしフォールバック(backboneは維持、
Cβのみ理想化配置)が使われています。変異結果は新しいstructureとして
保存され、元の`1a8o`は変更されません。

## 10. model(欠損残基検出・簡易ループ補完・ホモロジーモデリング)

```console
$ prot model gaps 1a8o
No gaps detected.
```

1A8Oには数値上の欠番がないので検出なし、というのは妥当な結果です。
欠番がある構造での`prot model loop`の挙動は、英語版チュートリアルの
Appendixで詳しく扱っています(直線補間による「粗いプレースホルダー」
であることを明記しています)。ホモロジーモデリングはMODELLERの外部
連携のみです:

```console
$ prot model homology --alignment a.pir --template t --target t2 \
    --template-dir . --output-dir out/
Error: The `modeller` Python package is not installed. Homology modeling
has no dependency-free fallback -- install MODELLER (license required
from https://salilab.org/modeller/) to use this command.
```

## 11. fix(PDBFixer連携での修復・正規化)

PDBFixerはpipで無料インストールできるので(`pip install -e ".[fix]"`)、
Scwrl4やMODELLERとは違って「常時使える」前提で組み込んであります。

```console
$ prot fix report 1a8o
p_bc8a53a4  (1a8o)
  missing residues: none found (PDBFixer needs SEQRES for this -- try
  `prot model gaps` too, which works from numbering alone)
  incomplete residues: none found
  nonstandard residues:
    A/MSE151 -> MET
    A/MSE185 -> MET
    A/MSE214 -> MET
    A/MSE215 -> MET

$ prot fix apply 1a8o --remove-heterogens water --name 1a8o_fixed
  residues completed (missing atoms added): 4
  nonstandard residues replaced: 4
  heterogens removed
Saved as '1a8o_fixed' (p_db5c7aa8)
```

セレノメチオニン(MSE)4残基がすべてMET(メチオニン)へ正規化され、
waterはそのまま保持されています。

**`prot model gaps`との役割分担**: PDBFixerの欠損残基検出はSEQRES
レコードが無いと機能しません。実際にこのプロジェクト自身の合成
フィクスチャ(`gapped.pdb`、番号が2→6に飛んでいて3残基欠損)で試すと、
PDBFixerは何も検出できませんが(SEQRESが無いため)、`prot model gaps`
は残基番号だけを見て正しく3残基欠損を検出します。**欠損検出は
`prot model gaps`、欠損原子の補完・非標準残基の正規化・heterogen除去・
水素付加は`prot fix`**、という使い分けです。

## 12. compare(構造比較)

```console
$ prot compare rmsd 1a8o 1a8o_A200THRALA
RMSD (fit, 70 common CA atoms): 0.000 A

$ prot compare secondary 1a8o 1a8o_A200THRALA
Secondary structure similarity: 0.89  (70 common residues)

$ prot compare contact 1a8o 1a8o_A200THRALA
Contact similarity (Jaccard): 1.00  (284/284 contacts shared)
```

`cb_only`変異はbackboneを一切変えないので、RMSDがきっちり0になるのは
妥当な健全性チェックです。secondary structure類似度が完全な1.0でない
のは、phi/psi分類器が側鎖変化による局所パッキングのわずかな違いにも
敏感に反応しているためです。`prot compare tmscore`/`pocket`/`ligand`
も同じ枠組みで使えます。

## 13. cluster(アンサンブルクラスタリング)

```console
$ prot cluster ensemble 1a8o 1a8o_A200THRALA --threshold 1.0
1 cluster(s) (method=greedy)
  #1: representative=1a8o + 1a8o_A200THRALA
```

RMSD 0の2構造は正しく同じクラスタにまとまります(このチュートリアルを
書く過程で、実はこの一致が壊れていたバグを実際に見つけて修正しました
— 詳細は英語版チュートリアルのAppendix参照)。`prot cluster models`は
1ファイル内の複数MODELをクラスタリングします。

## 14. annotate(アノテーション)

```console
$ prot annotate metadata 1a8o
method: x-ray diffraction
resolution: 1.7
deposition date: 1998-03-27
```

1A8Oには金属イオンが無いので、金属結合部位の検出は合成データ
(`examples/illustrative/metal_site.pdb`: Zn²⁺をHis/Cys/Aspが2.1Åで
配位)で実演します:

```console
$ prot import examples/illustrative/metal_site.pdb --name metalsite
$ prot annotate metal-sites metalsite
A/ZN100 (ZN): A/ASP3 (2.33 A), A/CYS2 (2.33 A), A/HIS1 (2.33 A)
```

`uniprot`/`pfam`は外部REST APIへの生きた接続が必要です。このチュート
リアルを書いたサンドボックス環境からは到達できず、エラーハンドリング
の実例として good demonstration になっています:

```console
$ prot annotate uniprot P12497
Error: HTTP 403 fetching https://rest.uniprot.org/uniprotkb/P12497.json: Forbidden
```

## 15. map(ビューア向けスクリプト生成)

```console
$ prot map mutation 1a8o mutation.pml --residue A/200 --tool pymol
Saved mutation.pml (1 residue(s), tool=pymol)
```

```pymol
# PyMOL coloring script (1 group(s))
color gray80, all
select grp_1, 1a8o and ((chain A and resi 200))
color red, grp_1
# grp_1 = mutations
```

`prot map pocket`/`domain`/`conservation`も同じ枠組みでPyMOL/ChimeraX/
VMD向けのスクリプトを生成します。

## 16. search(Foldseekによる構造類似性検索)

Foldseekはコンパイル済みバイナリ配布でPyPIには無いため、Scwrl4/
MODELLER/TMalignと同じ「外部ツールのみ・フォールバックなし」です。

```console
$ prot search foldseek 1a8o --against-project
Error: foldseek executable not found on PATH. Install Foldseek
(https://github.com/steineggerlab/foldseek#installation -- conda/homebrew/
static binary, not on PyPI) to search structural databases -- there is no
dependency-free fallback for large-scale structural similarity search.
```

`--against-project`はプロジェクト内の他構造をまとめて検索対象にできる
便利オプションです。`--target-db`/`--target-dir`で外部DBやディレクトリ
も指定できます。

## 17. predict / view(外部ツールのみのコマンド)

`prot predict`と`prot view`も、依存なし代替が原理的に成立しない領域
(構造予測モデル・GUIビューア起動)なので外部ツールのみです。

```console
$ prot predict colabfold "MDIRQ...ACQG" --output-dir out/
Error: colabfold_batch not found on PATH. Install ColabFold
(https://github.com/sokrypton/ColabFold) to predict structures locally --
there is no dependency-free fallback for structure prediction.

$ prot view 1a8o
Error: pymol executable not found on PATH. Install it from
https://pymol.org/ -- there is no dependency-free substitute for viewing
a structure.
```

## 18. valid(構造検証)

```console
$ prot valid clashes 1a8o
No clashes found.

$ prot valid geometry 1a8o
No bond geometry outliers found.
```

高分解能(1.7Å)の実データなので、原子重なり(clash)も結合幾何の逸脱も
ゼロ、という健全な結果です。どちらも依存なしで、vdW半径の重なりや
標準的な理想結合長・結合角(教科書的な値であり統計的較正値ではない)
という「誠実に検証可能な事実」だけを見ています。

Ramachandran outlier判定やrotamer outlier判定のような、統計的に較正
された参照データが必要な本格検証は、精度に自信が持てる保証がないため
自前実装せず、MolProbityへの外部連携のみにしています:

```console
$ prot valid molprobity 1a8o
Error: No MolProbity installation found on PATH (tried: phenix.molprobity,
mmtbx.molprobity, molprobity.molprobity). MolProbity ships as part of
Phenix (https://phenix-online.org/) or as a standalone build
(https://github.com/rlabduke/MolProbity) -- there is no dependency-free
substitute for its calibrated Ramachandran/rotamer/clashscore analysis;
see `prot valid clashes`/`prot valid geometry` for what this package can
check without it.
```

なお実装の途中、1A8Oの実データで**ジスルフィド結合(CYS198-CYS218)が
誤ってclash判定される**というバグを見つけて修正しました。共有結合を
単純な原子重なりと誤認していたのが原因で、ジスルフィド結合ペアを
除外するよう直しています。

## 19. replay(ワークフローの再実行)

プロジェクトに変更を加えるコマンドはすべて`.proteinexplorer/log.json`
に記録されます。`prot replay`はそのログを最初から再実行します。

```console
$ prot replay --dry-run
  [1] plan: import examples/1a8o/1A8O.pdb --name 1a8o
  [2] plan: mutate 1a8o --chain A --resid 200 --to ALA
  [3] plan: fix apply 1a8o --remove-heterogens water --name 1a8o_fixed
  [4] plan: import examples/illustrative/metal_site.pdb --name metalsite

$ prot replay
Backed up previous project state to .proteinexplorer_prereplay_20260814T050740
  [1] ok: import examples/1a8o/1A8O.pdb --name 1a8o
  [2] ok: mutate 1a8o --chain A --resid 200 --to ALA
  [3] ok: fix apply 1a8o --remove-heterogens water --name 1a8o_fixed
  [4] ok: import examples/illustrative/metal_site.pdb --name metalsite
4 step(s), 0 failed
```

structure IDはimportのたびに新しく生成されるので、後続コマンドの
argvに旧IDがリテラルで含まれていても、name照合で自動的に新IDへ
書き換えられます。ログを手で編集する必要はありません。

## 20. assembly(gemmi併用によるbiological assembly生成)

ここまでの全コマンドはBio.PDBというライブラリを土台にしていますが、
Bio.PDBには「結晶構造ファイルに書かれた対称操作(PDBのREMARK 350や
mmCIFの`_pdbx_struct_assembly_gen`)を展開して、実際の生物学的な複合体
(biological assembly)を組み立てる」機能がありません。非対称単位
(asymmetric unit、ファイルに実際に座標が書かれている範囲)しか扱えない
のです。これはプロジェクト最初期の仕様レビューで指摘したまま積み残し
になっていた課題で、gemmi(`pip install -e ".[assembly]"`、無料でpip
配布)を併用することで解決しました。

役割は徹底的に絞ってあります。gemmiは「assemblyを展開して書き出す」
という、この一点だけに使い、書き出した結果は普通のPDB/mmCIFファイル
として、これまで通りBio.PDBベースの`prot`の全コマンドでそのまま扱え
ます。

```console
$ prot assembly list 1a8o
  #1  DIMERIC  chains=A  operators=2
```

1A8Oのファイルには実は「これは二量体である」という情報が
(`DIMERIC`として)明記されています。展開してみます:

```console
$ prot assembly generate 1a8o
Assembly '1': ['A'] -> ['A1', 'A2']
Saved as '1a8o_assembly' (p_1ab81249)

$ prot info 1a8o_assembly
p_1ab81249  (1a8o_assembly)
  models: 1   chains: 2
  residues: 316   atoms: 1288
    protein: 140
    water: 176
```

非対称単位(chain A、644原子)から、対称操作を適用した実際の二量体
(chain A1+A2、1288原子)が正しく生成されました。生成された構造は
新しいstructureとしてプロジェクトに保存され、他の全コマンド
(`geometry`/`contact`/`compare`など)でそのまま使えます。

---

## まとめ: 三層構成の一覧

このツールキット全体を通じて、コマンドごとに以下の三層構成のどれかを
一貫して採用しています。

| パターン | 該当コマンド |
|---|---|
| **外部ツール優先+依存なしフォールバック** | `secondary`(DSSP/geometric)、`mutate`(Scwrl4/cb_only) |
| **常時利用可能なoptional extra(ライセンスフリー)** | `fix`(PDBFixer)、`assembly`(gemmi)、`cluster --method hierarchical`(scipy)、`plot`(matplotlib) |
| **外部ツールのみ・フォールバックなし** | `search`(Foldseek)、`predict`(ColabFold/AlphaFold)、`model homology`(MODELLER)、`view`(PyMOL/ChimeraX/VMD)、`valid molprobity`(MolProbity) |
| **依存なしの誠実な近似(検証可能な事実のみ)** | `pocket`(LIGSITE風グリッド探索)、`valid clashes`/`geometry`(vdW重なり・理想結合幾何) |

どの層を選ぶかは「そのツールがライセンスフリーでpip配布されているか」
「統計的に較正された参照データが必要か、それとも幾何学的に検証可能な
事実で足りるか」で判断してきました。精度に自信が持てないものは
無理に自前実装せず、正直に外部ツール任せにする、という方針を貫いて
います。
