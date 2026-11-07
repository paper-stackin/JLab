# Projected vs recreated \(\pi^+\) energy

Compares two ways of assigning the missing \(\pi^+\) energy in the CLAS12 exclusive channel

\[
e p \to e' p' \pi^- (\pi^+)
\]

and fills \(E_{\mathrm{rec}} - E_{\mathrm{proj}}\) in every \((W, Q^2)\) bin.

## Physics

Initial state (lab):

- electron \(e_i = (0,\,0,\,6.535,\,6.535)\) GeV
- proton \(p_i = (0,\,0,\,0,\,m_p)\)

The missing \(\pi^+\) four-momentum is reconstructed from four-momentum conservation:

\[
p_{\pi^+}^{\mathrm{miss}} = e_i + p_i - e_f - p_f - \pi^-
\]

| Quantity | Definition |
| --- | --- |
| Recreated energy | \(E_{\mathrm{rec}} = p_{\pi^+}^{\mathrm{miss}\,0}\) |
| Projected energy | \(E_{\mathrm{proj}} = \sqrt{\mathbf{p}_{\mathrm{miss}}^2 + m_{\pi}^2}\) (on-shell \(\pi^+\)) |

The histogram variable is \(E_{\mathrm{rec}} - E_{\mathrm{proj}}\). It is zero if the missing four-momentum is that of an on-shell \(\pi^+\).

### Frames

`use_cm_frame = false` (default): lab.

`use_cm_frame = true`: hadronic CM, i.e. the rest frame of the virtual photon plus the target proton, \(q + p_i\):

```cpp
TLorentzVector q = e_i - e_f;
TLorentzVector missing_pip = e_i + p_i - p_f - pi_minus - e_f;
missing_pip.Boost(-(q + p_i).BoostVector());
```

\(E_{\mathrm{rec}}\) and \(E_{\mathrm{proj}}\) are then taken from the boosted four-vector.

\(W\) and \(Q^2\) are always computed in the lab:

\[
Q^2 = -(e_i - e_f)^2,\qquad W = (e_i - e_f + p_i)^2
\]

## Input

Tree `MMpiptree` in

```
/home/stepan/downloads/missing_pip_clas12.root
```

Branches used: `px_e`, `py_e`, `pz_e`, `px_p`, `py_p`, `pz_p`, `px_pim`, `py_pim`, `pz_pim`.

Particle masses: \(m_p = 0.93827\) GeV, \(m_{\pi} = 0.13957\) GeV, \(m_e = 0\).

## Binning

| Axis | Edges / range |
| --- | --- |
| \(Q^2\) | \(0.5,\ 0.7,\ 1.0,\ 1.4,\ 2.0,\ 3.0,\ 5.5\) GeV\(^2\) (6 bins) |
| \(W\) | \(1.4\)–\(3.0\) GeV, \(25\) MeV bins (64 bins) |
| \(\Delta E\) | \([-0.3,\ 0.5]\) GeV, 100 bins |

Histogram names: `dE_Q2_bin=<q>_W_bin=<w>` with 1-based indices.

## Run

From this directory:

```bash
# Lab frame
root -l pi_plus_projected_vs_recreated_hists.cpp
root -l pi_plus_projected_vs_recreated_pdfs.cpp

# Hadronic CM frame
root -l 'pi_plus_projected_vs_recreated_hists.cpp(true)'
root -l 'pi_plus_projected_vs_recreated_pdfs.cpp(true)'
```

The first macro loops the full tree (~70M events) and should be run before the PDF macro.

## Output

| File | Contents |
| --- | --- |
| `pi_plus_projected_vs_recreated_hists.root` | Lab \(\Delta E\) histograms |
| `pi_plus_projected_vs_recreated_hists_cm.root` | CM \(\Delta E\) histograms |
| `dE_Q2_bin_N.pdf` | One PDF per \(Q^2\) bin (lab), one page per \(W\) bin |
| `dE_cm_Q2_bin_N.pdf` | Same for the CM frame |

## Macros

- `pi_plus_projected_vs_recreated_hists.cpp` — fill and write histograms
- `pi_plus_projected_vs_recreated_pdfs.cpp` — draw the ROOT file into PDFs
