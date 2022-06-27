@echo off
set layout=universal

for %%x in (A, B, C) do (
  echo %%x%layout%
  if %%x==A echo AAA
)

if %layout%==universal (
  echo hi%layout%
) else (
  echo bye
)