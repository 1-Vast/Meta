$ErrorActionPreference = 'Stop'
$cfg = Get-Content -Raw 'D:\MetaSieve\LLM\llm.json' | ConvertFrom-Json
$prompt = Get-Content -Raw 'D:\MetaSieve\LLM\k3_framework_prompt.txt'
if ($prompt.Length -lt 1) { throw 'K3 framework prompt was not loaded' }
$body = @{
    model = 'Kimi-K3'
    messages = @(@{ role = 'user'; content = $prompt })
    temperature = 0.2
    max_tokens = 700
} | ConvertTo-Json -Depth 8
$response = Invoke-RestMethod -Uri "$($cfg.llm.base_url)/chat/completions" -Method Post `
    -Headers @{ Authorization = "Bearer $($cfg.llm.api_key)" } `
    -ContentType 'application/json' -Body ([Text.Encoding]::UTF8.GetBytes($body)) `
    -TimeoutSec 700
$message = $response.choices[0].message
"# K3 Framework Review`n`n## Final Answer`n`n$($message.content)`n`n## Reasoning Transcript`n`n$($message.reasoning_content)" |
    Set-Content -LiteralPath 'D:\MetaSieve\report\K3_FRAMEWORK_REVIEW_20260814.md' -Encoding UTF8
