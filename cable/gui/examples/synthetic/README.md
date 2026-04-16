# Synthetic Examples

このディレクトリは公開文献ケースではなく，`cable_solver` と pycable GUI の入力形式・基本動作を確認するための synthetic / smoke-test 用 JSON を置きます．

| File | Format | Purpose |
|---|---|---|
| `catenary_500m.json` | legacy single-line | 静的カテナリの回帰確認．GUI integration test から参照されます． |
| `bridge_cable_small.json` | legacy single-line | 軽量な橋梁ケーブル scale の静的確認．GUI integration test から参照されます． |
| `dynamic_heave.json` | per-cable dynamic | 端点加振，時刻歴，snapshot 出力，GUI playback の確認． |
| `mooring_3leg_catenary.json` | multi-line BEM-compatible | `mooring_*` 13 要素配列形式の確認． |

公開文献に基づくベンチマークは一つ上の `examples/` 直下にあるケースディレクトリを使ってください．
