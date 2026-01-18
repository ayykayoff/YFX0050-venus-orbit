# YFX0050 – Keha liikumine Veenuse ümber (RKF45)

TalTech aine: **YMX0050 Füsikaliste protsesside modeleerimine**  
Õppejõud: **Mihhail Klopov**  

Autor: **Aimar Kamm**  
Üliõpilaskood: **233726IACB**

---

## Ülevaade
Selles töös simuleerin keha liikumist Veenuse ümber 2D tasandil, kasutades Runge–Kutta–Fehlbergi meetodit **RKF45**.

Simuleerin kaks juhtumit:
- **V = V1** (esimene kosmiline kiirus kõrgusel h0 = 200 km)
- **V = 1.2 · V1**

Joonistan:
- trajektoori **y(x)**
- energiad ajas **Ekin(t), Epot(t), Etot(t)** (per 1 kg)

---

## Failid
- `venus_orbit_rkf45.f90` — Fortrani põhiprogramm (kutsub RKF45 integraatorit)
- `srkf45.for` — RKF45 integraator (Sandia Labs)
- `orbit_v1.dat` — väljund (V = V1)
- `orbit_v12.dat` — väljund (V = 1.2·V1)
- `plot_two_cases.m` — MATLAB skript graafikute jaoks

---

## Kompileerimine ja käivitamine (Fortran)
```bash
gfortran -std=legacy venus_orbit_rkf45.f90 srkf45.for -o venus
./venus
