from pathlib import Path
SRC = Path("D:/MetaSieve/tools/research/stageX_csc_signal/stageQ2d1e_spaninit_interaction_20260818/runner_e.py")
DST = Path("D:/MetaSieve/tools/research/stageX_csc_signal/stageQ2d1e_spanparam_diag_20260819/runner_diag.py")
s = SRC.read_text(encoding="utf-8")
swaps = [
('''class InterOnly(nn.Module):
    def __init__(self, d_p, d_l):
        super().__init__()
        self.A = nn.Linear(d_p, RANK)
        self.B = nn.Linear(d_l, RANK)
        self.inter_bias = nn.Parameter(torch.zeros(1))
        self.inter_scale = nn.Parameter(torch.ones(1))

    def forward(self, p, l):
        inter = self.inter_scale * ((self.A(p) * self.B(l)).sum(-1) + self.inter_bias)
        return {"yhat": inter, "inter": inter}''',
'''class InterOnly(nn.Module):
    """Span-parameterized: A = V_train @ G with V_train fixed (diag prereg)."""
    def __init__(self, d_p, d_l):
        super().__init__()
        self.proj = nn.Linear(d_p, SPAN_RANK, bias=False)  # fixed V^T
        self.G = nn.Linear(SPAN_RANK, RANK)
        self.B = nn.Linear(d_l, RANK)
        self.inter_bias = nn.Parameter(torch.zeros(1))
        self.inter_scale = nn.Parameter(torch.ones(1))

    def forward(self, p, l):
        inter = self.inter_scale * ((self.G(self.proj(p)) * self.B(l)).sum(-1) + self.inter_bias)
        return {"yhat": inter, "inter": inter}'''),
('''        self.enc_p = nn.Linear(d_p, HID)
        self.enc_l = nn.Linear(d_l, HID)
        self.p_head = nn.Linear(HID, 1)
        self.l_head = nn.Linear(HID, 1)
        self.p_b = nn.Parameter(torch.zeros(n_rows))
        self.l_b = nn.Parameter(torch.zeros(n_lig))
        self.mu = nn.Parameter(torch.zeros(1))
        self.A = nn.Linear(HID, RANK)
        self.B = nn.Linear(HID, RANK)''',
'''        self.proj = nn.Linear(d_p, SPAN_RANK, bias=False)  # fixed V^T
        self.enc_p2 = nn.Linear(SPAN_RANK, HID)
        self.enc_l = nn.Linear(d_l, HID)
        self.p_head = nn.Linear(HID, 1)
        self.l_head = nn.Linear(HID, 1)
        self.p_b = nn.Parameter(torch.zeros(n_rows))
        self.l_b = nn.Parameter(torch.zeros(n_lig))
        self.mu = nn.Parameter(torch.zeros(1))
        self.G = nn.Linear(SPAN_RANK, RANK)
        self.B = nn.Linear(HID, RANK)'''),
('''        ep = self.enc_p(p)
        el = self.enc_l(l)
        pm = self.p_head(ep).squeeze(-1) + self.p_b[rows]
        lm = self.l_head(el).squeeze(-1) + self.l_b[ligs]
        inter = self.inter_scale * (((ep @ self.A.weight.T) * (el @ self.B.weight.T)).sum(-1)
                                    + self.inter_bias)''',
'''        ep = self.enc_p2(self.proj(p))
        el = self.enc_l(l)
        pm = self.p_head(ep).squeeze(-1) + self.p_b[rows]
        lm = self.l_head(el).squeeze(-1) + self.l_b[ligs]
        inter = self.inter_scale * ((self.G(self.proj(p)) * self.B(el)).sum(-1)
                                    + self.inter_bias)'''),
('''    if hasattr(model, "A"):
        with torch.no_grad():
            model.A.weight.mul_(0.5)
            if P.shape[1] == Vsp.shape[0] and Vsp is not None and arm != "oracle_diagnostic":
                # frozen span-initialized protein map: A = G @ Vsp^T, G xavier
                g = torch.empty(RANK, Vsp.shape[1])
                nn.init.xavier_uniform_(g)
                model.A.weight.copy_(g.to(model.A.weight.device)
                                     @ torch.from_numpy(Vsp.T).float().to(model.A.weight.device))
    if hasattr(model, "enc_p") and Vsp is not None and P.shape[1] == Vsp.shape[0]:
        with torch.no_grad():
            proj = torch.from_numpy((Vsp @ Vsp.T).astype(np.float32)).to(model.enc_p.weight.device)
            model.enc_p.weight.copy_(proj)''',
'''    if hasattr(model, "proj") and Vsp is not None:
        with torch.no_grad():
            model.proj.weight.copy_(torch.from_numpy(Vsp.T.astype(np.float32)).to(model.proj.weight.device))
            model.proj.weight.requires_grad_(False)
    if hasattr(model, "G"):
        with torch.no_grad():
            model.G.weight.mul_(0.5)'''),
('''        if hasattr(model, "A"):
            pen = L2_PEN * (model.A.weight.square().sum()
                            + model.B.weight.square().sum())''',
'''        if hasattr(model, "G"):
            pen = L2_PEN * (model.G.weight.square().sum()
                            + model.B.weight.square().sum())'''),
('''        res[arm] = eval_arm(best_model, P, arm, t, level, splits, device, Lt_dev)
        print(mech, level, seed, arm,''',
'''        res[arm] = eval_arm(best_model, P, arm, t, level, splits, device, Lt_dev)
        if arm == "correct" and hasattr(best_model, "G") and Vsp is not None:
            res[arm]["span_energy"] = 1.0  # by construction: A = V_train @ G
        print(mech, level, seed, arm,'''),
('out = {"schema": "MetaSieve.StageQ2d1e.LADDER.v1",',
 'out = {"schema": "MetaSieve.StageQ2d1eSpanParam.LADDER.v1",'),
('json.dump(out, open(HERE / "Q2D1E_LADDER.json", "w"), indent=1)',
 'json.dump(out, open(HERE / "Q2D1E_SPANPARAM_LADDER.json", "w"), indent=1)'),
('''RANK = 4
HID = 32''',
'''RANK = 4
SPAN_RANK = 28  # train-row feature span rank (frozen in Q2d-1c forensics)
HID = 32'''),
]
missing = []
for a, b in swaps:
    if a not in s:
        missing.append(a[:60].replace("\n", " | "))
    else:
        s = s.replace(a, b, 1)
if missing:
    print("MISSING", len(missing))
    for m in missing:
        print(" -", m)
else:
    DST.write_text(s, encoding="utf-8")
    print("wrote", DST, len(s))
