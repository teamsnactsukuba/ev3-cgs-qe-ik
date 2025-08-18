# algn

sqrt(2) を代数的数として扱った計算を記録

## スクリプト

* cgs-F-input-algn.0.rr: 単項式順序 DegRevLex でCGSを計算するスクリプト。出力先はG0.dat。
* cgs-F-input-algn.2.rr: 単項式順序 Lex でCGSを計算するスクリプト。出力先はG2.dat。
* extract-cgs-G.rr: Nabeshima による CGS 計算プログラム kcgs_main の出力から断片とGroebner基底を抽出するスクリプト。
* extract-cgs-g-script.rr: extract-cgs-G.rr を用いて CGS から断片とGroebner基底を抽出するスクリプト。G2-reverse.dat を用いて逆順（
* zero-dimensional-test.rr: CGS の(断片,基底)のペアに収められている基底 G で生成されるイデアル \<G\> がゼロ次元かどうかを判定する。実際には G に属する多項式 g の頭単項式
 (leading monomial) が1変数であり、すべての変数が網羅されているかどうかを調べる。

## データファイル

* G0.dat: 単項式順序 DegRevLex で計算したCGS. (断片,基底)のペアは計算順に収められている。成分のペアは最後に計算されたものが最初に入っているので注意。
* G2.dat: 辞書式順序 Lex で計算したCGS. (断片,基底)のペアは計算順に収められている。成分のペアは最後に計算されたものが最初に入っているので注意。
* G2-reverse.dat: G2.dat の(断片,基底)のペアの並び方を逆にしたもの。成分のペアのは計算した順番に収められている。