set layout=universal
chcp 65001
for %%G in (A, B, C) do (
  echo %%G%layout%
  if %%G==A echo AAA
)

if %layout%==universal (
  echo hi%layout%
) else (
  echo bye
)

echo "Δ"