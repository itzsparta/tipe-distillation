# TIPE : optimisation thermodynamique d'un alambic

Simulations Python de mon TIPE de CPGE TSI (2025–2026). Elles modélisent une
distillation discontinue éthanol-eau pour répondre à une question : comment
régler la puissance de chauffe pour obtenir un distillat plus pur sans y passer
des heures ?

## Modèle

- **Équilibre liquide-vapeur** : courbe expérimentale éthanol-eau (fraction
  molaire liquide x, vapeur y, température d'ébullition), interpolée
  linéairement.
- **Bilan matière pas à pas** (dt = 5 s) : la puissance de chauffe, divisée par
  l'enthalpie de vaporisation du mélange, donne les moles vaporisées. La vapeur
  a la composition d'équilibre y(x) et appauvrit la cuve en éthanol.
- **Conversion** des fractions molaires en degré alcoolique volumique (% ABV) à
  partir des masses molaires et des masses volumiques.
- **Coupes** : têtes (premiers 200 mL), cœur, queues (sous un seuil de degré).

## Scripts

| Script | Question |
|---|---|
| `simulation_coupes.py` | 15 L de moût à 12 %, chauffe en escalier de 800 à 2000 W. Courbes du degré du distillat (têtes, cœur, queues), du débit et de la puissance. |
| `distillations_successives.py` | Gain de pureté de trois passes successives (2500, 1500 puis 800 W), chacune repartant du distillat de la précédente. |
| `optimisation_puissance.py` | Compromis puissance / pureté / durée, de 800 à 2400 W, imprimé sous forme de tableau. |

## Limites

- Un seul étage d'équilibre, sans reflux : l'alambic est traité comme un
  bouilleur simple.
- La perte d'efficacité au-delà de 800 W dans `optimisation_puissance.py` est
  une hypothèse linéaire, pas une mesure.
- Pas de pertes thermiques : toute la puissance sert à vaporiser.

## Lancer

```bash
pip install -r requirements.txt
python simulation_coupes.py
```
