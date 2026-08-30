import { describe, expect, it } from 'vitest';

import { buildOrbitalPresentation } from '../../src/app/orbitalPresentation';
import {
  hydrogenSchrodingerEnergyElectronVolts,
  HYDROGEN_REDUCED_MASS_RATIO,
} from '../../src/science/hydrogen/energy';
import {
  hydrogenExpectedRadiusBohr,
  hydrogenNodeCounts,
} from '../../src/science/hydrogen/observables';
import type { RealOrbitalName } from '../../src/science/hydrogen/realOrbitals';
import {
  bohrToMeters,
  bohrToPicometers,
  metersToNanometers,
} from '../../src/science/units/atomicUnits';
import type { OrbitalSamplingState } from '../../src/sampling/contracts';

describe('présentation scientifique d’un état orbital', () => {
  it('présente un état complexe sans l’assimiler à une orbitale réelle', () => {
    const state = { basis: 'complex', n: 3, l: 2, m: 2 } as const;
    const presentation = buildOrbitalPresentation(state);
    const expectedRadiusBohr = hydrogenExpectedRadiusBohr(3, 2);

    expect(presentation).toEqual({
      angularDegree: 2,
      basisLabel: 'Base complexe |n,l,m⟩',
      energy: {
        electronVolts: hydrogenSchrodingerEnergyElectronVolts(3),
        modelLabel: '¹H non relativiste — masse réduite électron-proton',
        reducedMassRatio: HYDROGEN_REDUCED_MASS_RATIO,
        unit: 'eV',
      },
      expectedRadius: {
        bohr: expectedRadiusBohr,
        nanometers: metersToNanometers(bohrToMeters(expectedRadiusBohr)),
        picometers: bohrToPicometers(expectedRadiusBohr),
        unitLabels: {
          bohr: 'a₀',
          nanometers: 'nm',
          picometers: 'pm',
        },
      },
      nodes: hydrogenNodeCounts(3, 2),
      notation: '3d (m = +2)',
      realCombinationLabel: null,
    });
  });

  it('dérive le degré et la combinaison ±m d’une orbitale réelle', () => {
    const presentation = buildOrbitalPresentation({ basis: 'real', n: 3, orbital: 'd_xy' });

    expect(presentation.angularDegree).toBe(2);
    expect(presentation.notation).toBe('3d_xy');
    expect(presentation.basisLabel).toBe('Base réelle normalisée');
    expect(presentation.realCombinationLabel).toBe(
      'Combinaison réelle normalisée des états m = ±2',
    );
    expect(presentation.nodes).toEqual({ radialNodes: 0, angularNodes: 2, totalNodes: 2 });
  });

  it('décrit séparément une harmonique zonale réelle m=0', () => {
    const presentation = buildOrbitalPresentation({ basis: 'real', n: 2, orbital: 'p_z' });

    expect(presentation.notation).toBe('2p_z');
    expect(presentation.realCombinationLabel).toBe('État réel m = 0 (harmonique zonale)');
  });

  it('conserve une notation non ambiguë au-delà des lettres prises en charge', () => {
    const presentation = buildOrbitalPresentation({ basis: 'complex', n: 9, l: 8, m: -7 });

    expect(presentation.notation).toBe('n=9, l=8 (m = -7)');
    expect(presentation.nodes).toEqual({ radialNodes: 0, angularNodes: 8, totalNodes: 8 });
  });

  it('rejette les états complexes et réels invalides', () => {
    const invalidComplex = { basis: 'complex', n: 2, l: 2, m: 0 } as OrbitalSamplingState;
    const invalidRealName = {
      basis: 'real',
      n: 3,
      orbital: 'f_xyz' as RealOrbitalName,
    } as const;
    const invalidRealPrincipalNumber = { basis: 'real', n: 2, orbital: 'd_xy' } as const;

    expect(() => buildOrbitalPresentation(invalidComplex)).toThrow(RangeError);
    expect(() => buildOrbitalPresentation(invalidRealName)).toThrow(RangeError);
    expect(() => buildOrbitalPresentation(invalidRealPrincipalNumber)).toThrow(RangeError);
  });
});
