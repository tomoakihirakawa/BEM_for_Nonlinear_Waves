# Placeholder cables — Fred Hartman Bridge

`C_short.json`, `C_mid.json`, `C_long.json` は Fred Hartman Bridge の
実測端点座標ではない**プレースホルダ**です．代表的な短・中・長スパン
1 本ずつを便宜的に配置しただけで，FHWA-HRT-05-083 の特定ケーブルに
対応しません．

実データの検証は親ディレクトリの [published_table15/](../published_table15/)
(Table 15, 13S-24S の 12 本) を使ってください．そちらは published
mass/T/L/f の内部整合と solver equilibrium の対比ができます．

## 用途

- GUI / solver の **smoke test** (3 本で dynamic mode を軽く走らせる)
- 新機能の回帰確認

published validation には使わないでください．
