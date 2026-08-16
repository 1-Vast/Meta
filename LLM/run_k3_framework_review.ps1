$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.Net.Http
$cfg = Get-Content -Raw 'D:\MetaSieve\LLM\llm.json' | ConvertFrom-Json
$prompt = Get-Content -Raw 'D:\MetaSieve\LLM\k3_framework_prompt.txt'
if ($prompt.Length -lt 500) { throw 'K3 framework prompt was not loaded' }
$body = @{
    model = 'Kimi-K3'
    messages = @(@{ role = 'user'; content = $prompt })
    temperature = 0.2
    max_tokens = 4500
    stream = $true
} | ConvertTo-Json -Depth 8 -Compress
$client = [System.Net.Http.HttpClient]::new()
$client.Timeout = [TimeSpan]::FromMinutes(12)
$client.DefaultRequestHeaders.Authorization = [System.Net.Http.Headers.AuthenticationHeaderValue]::new('Bearer', $cfg.llm.api_key)
$request = [System.Net.Http.HttpRequestMessage]::new([System.Net.Http.HttpMethod]::Post, "$($cfg.llm.base_url)/chat/completions")
$request.Content = [System.Net.Http.StringContent]::new($body, [Text.Encoding]::UTF8, 'application/json')
$response = $client.SendAsync($request, [System.Net.Http.HttpCompletionOption]::ResponseHeadersRead).GetAwaiter().GetResult()
$response.EnsureSuccessStatusCode()
$reader = [IO.StreamReader]::new($response.Content.ReadAsStreamAsync().GetAwaiter().GetResult())
$reasoning = [Text.StringBuilder]::new()
$content = [Text.StringBuilder]::new()
while (($line = $reader.ReadLine()) -ne $null) {
    if (-not $line.StartsWith('data: ') -or $line -eq 'data: [DONE]') { continue }
    $chunk = $line.Substring(6) | ConvertFrom-Json
    $delta = $chunk.choices[0].delta
    if ($delta.reasoning_content) { [void]$reasoning.Append($delta.reasoning_content) }
    if ($delta.content) { [void]$content.Append($delta.content) }
}
"# K3 Framework Review`n`n## Final Answer`n`n$content`n`n## Reasoning Transcript`n`n$reasoning" |
    Set-Content -LiteralPath 'D:\MetaSieve\report\K3_FRAMEWORK_REVIEW_20260814.md' -Encoding UTF8
