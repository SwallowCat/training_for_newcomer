# 配列情報処理

ヒト21番染色体の先頭コンティグのDNA配列ファイル `data/NT_113952.1.fasta` ファイルを用いて以下の課題を実行せよ。

## 1-1 塩基数カウント

与えられたDNA配列に含まれるA/T/G/Cの各塩基数を出力せよ。

## 1-2 逆相補鎖の作成

与えられたDNA配列の逆相補鎖を出力せよ。

## 1-3 GC含量推移プロット

ウィンドウ幅 $w$、ステップ幅 $s$ で、ウィンドウ内のGCの割合（百分率）を連続的に出力するプログラムを作成し、 $w=1000, s=300$ で得られる値を出力せよ。また、その結果をmatplotlibを用いて図示せよ

<img src="https://user-images.githubusercontent.com/6902135/229272441-b245e68d-7679-4351-b46a-8c8a8000de8a.png" style="width:50%">

## 1-4 塩基配列検索

引数で与えられた部分配列（クエリ）を検索し、その部分配列が出現する場所をすべて列挙せよ。ただし、順方向だけではなく逆相補鎖上も検索すること。

列挙の際には、その部分配列の開始点を1-originで出力し、順方向の場合は `F`、逆相補鎖の場合は `R` を接頭辞としてつけること。逆相補鎖における「開始点」は塩基配列上最も右側になる。

例
```
配列　: atgccgt
クエリ: cg
出力　: F5, R6
```

## 1-5 アミノ酸配列への翻訳

塩基配列をアミノ酸配列（1文字表記の列）に翻訳せよ。
出力されるアミノ酸配列は、翻訳開始位置（これは必ずメチオニン `M` になる）からStopコドン `_` までで1行の出力を構成し、次の翻訳開始位置からは次の行に記載すること。

なお、翻訳の際には1つの配列に対して3つの読み枠（リーディングフレーム）が存在し、さらに塩基配列は逆相補鎖があるので、合計6つの読み枠に対して上記の処理を行うこと。