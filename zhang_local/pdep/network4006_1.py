species(
    label = '[CH]C(=C)C[C]([O])[O](18595)',
    structure = SMILES('[CH]C(=C)C[C]([O])[O]'),
    E0 = (611.706,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,360,370,350,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31777,0.0684171,-0.000101379,9.96189e-08,-3.90509e-11,73658.6,29.1468], Tmin=(100,'K'), Tmax=(814.487,'K')), NASAPolynomial(coeffs=[1.58731,0.0464918,-2.30589e-05,4.45808e-09,-3.09959e-13,74298,32.0967], Tmin=(814.487,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(611.706,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Cs_P) + radical(AllylJ2_triplet) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = 'C=C([O])[O](1172)',
    structure = SMILES('C=C([O])[O]'),
    E0 = (-40.8548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,350,440,435,1725,219.703,219.72],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3790.78,'J/mol'), sigma=(6.02099,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.11 K, Pc=39.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.85742,0.0201906,-6.55878e-06,-8.62627e-09,5.39429e-12,-4868.13,11.9276], Tmin=(100,'K'), Tmax=(974.512,'K')), NASAPolynomial(coeffs=[9.20068,0.00568109,-1.96826e-06,3.71357e-10,-2.78204e-14,-6651.8,-21.3196], Tmin=(974.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.8548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'C#C[CH2](17441)',
    structure = SMILES('C#C[CH2]'),
    E0 = (328.481,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2175,525,1131.03,1132.16,1135.9],'cm^-1')),
        HinderedRotor(inertia=(0.154206,'amu*angstrom^2'), symmetry=1, barrier=(3.5455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2095.25,'J/mol'), sigma=(4.76,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32026,0.0108736,8.62061e-06,-1.82973e-08,7.68649e-12,39535.3,8.27851], Tmin=(100,'K'), Tmax=(960.555,'K')), NASAPolynomial(coeffs=[6.38511,0.00814486,-2.78734e-06,4.95348e-10,-3.50148e-14,38483.6,-8.79383], Tmin=(960.555,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Propargyl)"""),
)

species(
    label = '[CH]C1([CH2])CC1([O])[O](20876)',
    structure = SMILES('[CH]C1([CH2])CC1([O])[O]'),
    E0 = (705.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25877,0.0525825,-3.81442e-05,1.32203e-08,-1.81567e-12,85013.3,25.7343], Tmin=(100,'K'), Tmax=(1709.89,'K')), NASAPolynomial(coeffs=[15.9101,0.0183079,-8.07651e-06,1.49714e-09,-1.01627e-13,80002.9,-52.8092], Tmin=(1709.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(705.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsCs) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)(O)OJ) + radical(Neopentyl) + radical(CCJ2_triplet) + radical(CC(C)(O)OJ)"""),
)

species(
    label = '[CH]C1([CH2])C[C]([O])O1(20824)',
    structure = SMILES('[CH]C1([CH2])C[C]([O])O1'),
    E0 = (714.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.769298,0.0678644,-8.01395e-05,5.08475e-08,-1.24443e-11,86106.4,26.437], Tmin=(100,'K'), Tmax=(1124.3,'K')), NASAPolynomial(coeffs=[12.6914,0.0184066,-4.75968e-06,5.79186e-10,-2.78643e-14,83870.7,-30.4976], Tmin=(1124.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(714.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCOJ) + radical(CJC(C)OC) + radical(CCJ2_triplet) + radical(Cs_P)"""),
)

species(
    label = 'H(8)',
    structure = SMILES('[H]'),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,25474.2,-0.444973], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]C(=C)[CH]C([O])=O(20877)',
    structure = SMILES('[CH]C(=C)[CH]C([O])=O'),
    E0 = (315.426,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58156,0.0495455,-2.41053e-05,2.39998e-09,5.60463e-13,38027,21.6248], Tmin=(100,'K'), Tmax=(1623.35,'K')), NASAPolynomial(coeffs=[16.2744,0.0255357,-1.31872e-05,2.54344e-09,-1.74232e-13,31650,-61.3268], Tmin=(1623.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(315.426,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(AllylJ2_triplet) + radical(C=CCJCO)"""),
)

species(
    label = '[CH][C]=C(18825)',
    structure = SMILES('[CH][C]=C'),
    E0 = (614.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,228.264,228.889,229.07],'cm^-1')),
        HinderedRotor(inertia=(1.35219,'amu*angstrom^2'), symmetry=1, barrier=(50.6528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27541,0.0127954,9.49515e-06,-1.56026e-08,5.42938e-12,73954,11.3502], Tmin=(100,'K'), Tmax=(1063.31,'K')), NASAPolynomial(coeffs=[4.18965,0.0168435,-6.77763e-06,1.22218e-09,-8.33556e-14,73336.3,4.89309], Tmin=(1063.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)[CH]C([O])[O](20878)',
    structure = SMILES('[CH]C(=C)[CH]C([O])[O]'),
    E0 = (523.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10639,0.0711035,-9.67834e-05,8.88526e-08,-3.41253e-11,63044.6,26.4181], Tmin=(100,'K'), Tmax=(778.606,'K')), NASAPolynomial(coeffs=[3.69307,0.0451744,-2.24783e-05,4.37933e-09,-3.07122e-13,63024.9,17.0466], Tmin=(778.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(523.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(CCOJ) + radical(CCOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]C(=C)[CH][C]([O])O(20879)',
    structure = SMILES('[CH]C(=C)[CH][C]([O])O'),
    E0 = (502.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.972377,0.0720425,-8.45092e-05,5.98356e-08,-1.81224e-11,60591.1,27.4375], Tmin=(100,'K'), Tmax=(788.209,'K')), NASAPolynomial(coeffs=[7.8732,0.0370226,-1.78653e-05,3.46883e-09,-2.445e-13,59503.2,-4.21268], Tmin=(788.209,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(AllylJ2_triplet) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[CH]C(=[CH])CC([O])[O](20880)',
    structure = SMILES('[CH]C(=[CH])CC([O])[O]'),
    E0 = (653.556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,350,440,435,1725,455.441,455.452,455.453,455.454,455.456,455.463],'cm^-1')),
        HinderedRotor(inertia=(0.000812587,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000812648,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.355424,'amu*angstrom^2'), symmetry=1, barrier=(52.3211,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1999,0.0712805,-0.000107796,1.04917e-07,-4.05337e-11,78696,28.5923], Tmin=(100,'K'), Tmax=(819.256,'K')), NASAPolynomial(coeffs=[2.11198,0.0458197,-2.27159e-05,4.38361e-09,-3.04178e-13,79251.5,28.6765], Tmin=(819.256,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(AllylJ2_triplet) + radical(CCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=[CH])C[C]([O])O(20881)',
    structure = SMILES('[CH]C(=[CH])C[C]([O])O'),
    E0 = (633.097,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,360,370,350,3615,1277.5,1000,350,440,435,1725,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.912865,0.0742637,-0.000103952,8.90347e-08,-3.13228e-11,76249,30.1461], Tmin=(100,'K'), Tmax=(800.033,'K')), NASAPolynomial(coeffs=[6.55637,0.0371828,-1.78077e-05,3.40064e-09,-2.35371e-13,75629.7,5.9516], Tmin=(800.033,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(633.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Cs_P) + radical(CCOJ) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C]([O])[O](2304)',
    structure = SMILES('[CH2][C]([O])[O]'),
    E0 = (411.444,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3000,3100,440,815,1455,1000,1740.05,1740.26],'cm^-1')),
        HinderedRotor(inertia=(0.00349141,'amu*angstrom^2'), symmetry=1, barrier=(7.50228,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01691,0.0315826,-7.13533e-05,8.28598e-08,-3.37644e-11,49511,16.8536], Tmin=(100,'K'), Tmax=(857.006,'K')), NASAPolynomial(coeffs=[-0.765291,0.0236786,-1.27872e-05,2.50409e-09,-1.7282e-13,51097.8,39.9925], Tmin=(857.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CJCO) + radical(CCOJ) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[CH][C]([CH2])C=C([O])[O](20882)',
    structure = SMILES('[CH][C]([CH2])C=C([O])[O]'),
    E0 = (654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,977.413,977.481,977.507,977.568],'cm^-1')),
        HinderedRotor(inertia=(0.00499504,'amu*angstrom^2'), symmetry=1, barrier=(3.38755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147329,'amu*angstrom^2'), symmetry=1, barrier=(3.38737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0918885,'amu*angstrom^2'), symmetry=1, barrier=(62.3416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09594,0.0598679,-6.15959e-05,3.13741e-08,-5.7943e-12,78766.2,26.6994], Tmin=(100,'K'), Tmax=(963.395,'K')), NASAPolynomial(coeffs=[13.5198,0.0160518,-5.46867e-06,9.03389e-10,-5.89556e-14,76011.9,-34.6461], Tmin=(963.395,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(C=COJ) + radical(Allyl_T) + radical(C=COJ)"""),
)

species(
    label = '[CH]C(=[CH])C[C]([O])[O](20883)',
    structure = SMILES('[CH]C(=[CH])C[C]([O])[O]'),
    E0 = (858.802,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,360,370,350,412.494,412.51,412.51,412.539,412.542,412.563],'cm^-1')),
        HinderedRotor(inertia=(0.000990659,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000990658,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000990657,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2205,0.0722086,-0.000118897,1.17903e-07,-4.51986e-11,103379,29.2061], Tmin=(100,'K'), Tmax=(833.792,'K')), NASAPolynomial(coeffs=[2.0073,0.0433468,-2.18421e-05,4.21565e-09,-2.91439e-13,104120,30.7826], Tmin=(833.792,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(858.802,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cs_P) + radical(Cds_P) + radical(AllylJ2_triplet) + radical(CCOJ)"""),
)

species(
    label = '[CH][C]1CC([O])([O])C1(20884)',
    structure = SMILES('[CH][C]1CC([O])([O])C1'),
    E0 = (694.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64661,0.0387727,-1.77586e-05,2.78854e-09,-6.87267e-14,83512.9,22.1535], Tmin=(100,'K'), Tmax=(2243.6,'K')), NASAPolynomial(coeffs=[20.2994,0.0153249,-7.44705e-06,1.3187e-09,-8.25752e-14,73572.1,-81.7768], Tmin=(2243.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(694.032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)(O)OJ) + radical(Tertalkyl) + radical(CC(C)(O)OJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C]1CO[C]([O])C1(20885)',
    structure = SMILES('[CH][C]1CO[C]([O])C1'),
    E0 = (600.126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23673,0.028771,2.62734e-05,-5.41406e-08,2.35568e-11,72251.2,22.4115], Tmin=(100,'K'), Tmax=(913.467,'K')), NASAPolynomial(coeffs=[9.5156,0.0228746,-6.70131e-06,1.05701e-09,-7.00386e-14,69837.6,-17.9784], Tmin=(913.467,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.126,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(Cs_P) + radical(CCOJ) + radical(CCJ(C)CO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])[CH][C]([O])[O](20886)',
    structure = SMILES('[CH]C([CH2])[CH][C]([O])[O]'),
    E0 = (976.019,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,360,370,350,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,232.152,1238.09,1380.18,1447.36,1529.35,1904.57,2180.82],'cm^-1')),
        HinderedRotor(inertia=(0.0737854,'amu*angstrom^2'), symmetry=1, barrier=(2.83312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0737854,'amu*angstrom^2'), symmetry=1, barrier=(2.83312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0737854,'amu*angstrom^2'), symmetry=1, barrier=(2.83312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0737854,'amu*angstrom^2'), symmetry=1, barrier=(2.83312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 8,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26499,0.0693843,-0.00011228,1.05976e-07,-3.87441e-11,117478,34.2456], Tmin=(100,'K'), Tmax=(848.792,'K')), NASAPolynomial(coeffs=[4.08864,0.0358484,-1.72651e-05,3.2701e-09,-2.23403e-13,117727,25.3786], Tmin=(848.792,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(976.019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(CCJCO) + radical(Cs_P) + radical(Isobutyl) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH]C(=C)C([CH2])([O])[O](18614)',
    structure = SMILES('[CH]C(=C)C([CH2])([O])[O]'),
    E0 = (608.958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,180,1005.14,1005.14,1005.14,1005.14,1005.14,1005.14,1005.14,1005.14,1005.14,2307.69],'cm^-1')),
        HinderedRotor(inertia=(0.0577144,'amu*angstrom^2'), symmetry=1, barrier=(1.32697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0577144,'amu*angstrom^2'), symmetry=1, barrier=(1.32697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0577144,'amu*angstrom^2'), symmetry=1, barrier=(1.32697,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.407079,0.0800563,-9.73288e-05,5.98523e-08,-1.34083e-11,73369.3,28.3582], Tmin=(100,'K'), Tmax=(797.529,'K')), NASAPolynomial(coeffs=[13.9069,0.022518,-8.23768e-06,1.36862e-09,-8.75738e-14,70892.6,-35.744], Tmin=(797.529,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(608.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(C=CC(C)(O)OJ) + radical(C=CC(O)2CJ) + radical(C=CC(C)(O)OJ)"""),
)

species(
    label = '[CH]C(=C)CC([O])=O(18773)',
    structure = SMILES('[CH]C(=C)CC([O])=O'),
    E0 = (198.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84152,0.0451159,-1.47529e-05,-3.47571e-09,1.81116e-12,23954.6,22.7981], Tmin=(100,'K'), Tmax=(1564.81,'K')), NASAPolynomial(coeffs=[13.5237,0.0304234,-1.52104e-05,2.91444e-09,-1.99536e-13,18441.3,-44.7271], Tmin=(1564.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[O][C][O](2319)',
    structure = SMILES('[O][C][O]'),
    E0 = (504.804,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([786.798,787.166,787.343],'cm^-1')),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (44.0095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.36312,0.000897068,5.83704e-06,-5.0699e-09,1.05455e-12,60693,5.58445], Tmin=(100,'K'), Tmax=(1870.77,'K')), NASAPolynomial(coeffs=[6.74008,0.00434757,-3.77128e-06,7.92205e-10,-5.4643e-14,58310.5,-11.3625], Tmin=(1870.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsHH) + radical(CH2_triplet) + radical(OCOJ) + radical(OCOJ)"""),
)

species(
    label = '[CH]C([CH2])=C(18776)',
    structure = SMILES('[CH]C([CH2])=C'),
    E0 = (489.252,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,340.08,340.579,341.119],'cm^-1')),
        HinderedRotor(inertia=(0.611693,'amu*angstrom^2'), symmetry=1, barrier=(50.5163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.611621,'amu*angstrom^2'), symmetry=1, barrier=(50.5272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48527,0.0254697,1.30556e-05,-3.21059e-08,1.33815e-11,58905,14.7489], Tmin=(100,'K'), Tmax=(981.423,'K')), NASAPolynomial(coeffs=[8.34326,0.0203037,-7.64415e-06,1.37964e-09,-9.63641e-14,56854.1,-17.9932], Tmin=(981.423,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'CH2(T)(28)',
    structure = SMILES('[CH2]'),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154979,3.26298e-06,-2.40422e-09,5.69497e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.58,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50338e-15,46058.1,4.77808], Tmin=(1104.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C#CC[C]([O])[O](5453)',
    structure = SMILES('C#CC[C]([O])[O]'),
    E0 = (446.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,2175,525,360,370,350,180,1228.01,1229.39],'cm^-1')),
        HinderedRotor(inertia=(0.0960217,'amu*angstrom^2'), symmetry=1, barrier=(2.20773,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0976735,'amu*angstrom^2'), symmetry=1, barrier=(2.24571,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4078.22,'J/mol'), sigma=(6.60844,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=637.01 K, Pc=32.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07131,0.0527101,-9.72356e-05,9.92512e-08,-3.7486e-11,53722.6,22.4214], Tmin=(100,'K'), Tmax=(868.79,'K')), NASAPolynomial(coeffs=[1.58746,0.0299224,-1.47016e-05,2.77664e-09,-1.88028e-13,54750.8,30.1209], Tmin=(868.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH][C](C)C=C([O])[O](20887)',
    structure = SMILES('[CH][C](C)C=C([O])[O]'),
    E0 = (448.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924984,0.0598987,-5.46536e-05,2.54052e-08,-4.66249e-12,54109.5,24.8434], Tmin=(100,'K'), Tmax=(1323.35,'K')), NASAPolynomial(coeffs=[14.9951,0.0173699,-6.4479e-06,1.12055e-09,-7.47937e-14,50385.6,-46.9789], Tmin=(1323.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(448.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(CCJ2_triplet) + radical(Allyl_T) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2][C]([CH2])C=C([O])[O](18594)',
    structure = SMILES('[CH2][C]([CH2])C=C([O])[O]'),
    E0 = (410.867,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,478.01,478.097,478.121],'cm^-1')),
        HinderedRotor(inertia=(0.000737616,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000737894,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000737352,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23104,0.0576627,-5.12773e-05,1.92509e-08,-5.72106e-13,49518.8,27.6897], Tmin=(100,'K'), Tmax=(868.027,'K')), NASAPolynomial(coeffs=[12.3603,0.0188609,-5.79735e-06,8.88866e-10,-5.5313e-14,47116.4,-27.1361], Tmin=(868.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(Isobutyl) + radical(Allyl_T) + radical(C=COJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C(C)[CH]C([O])=O(20888)',
    structure = SMILES('[CH]=C(C)[CH]C([O])=O'),
    E0 = (191.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8168,0.0497349,-3.02609e-05,7.60094e-09,-7.16413e-13,23149.7,20.059], Tmin=(100,'K'), Tmax=(2517.84,'K')), NASAPolynomial(coeffs=[27.6438,0.00870595,-5.81892e-06,1.12949e-09,-7.38767e-14,10143.5,-128.391], Tmin=(2517.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CCJCO) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(=C)[CH]C([O])=O(18587)',
    structure = SMILES('[CH2]C(=C)[CH]C([O])=O'),
    E0 = (96.2408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42231,0.0491341,-2.3491e-05,2.97639e-11,1.66509e-12,11674,21.3578], Tmin=(100,'K'), Tmax=(1383.65,'K')), NASAPolynomial(coeffs=[15.0766,0.0244155,-1.26893e-05,2.53223e-09,-1.79556e-13,6483.04,-54.0539], Tmin=(1383.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.2408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + radical(Allyl_P) + radical(CCOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]=[C]CC[C]([O])[O](18582)',
    structure = SMILES('[CH]=[C]CC[C]([O])[O]'),
    E0 = (741.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3120,650,792.5,1650,360,370,350,180,180,1964.91,1965.18],'cm^-1')),
        HinderedRotor(inertia=(0.12083,'amu*angstrom^2'), symmetry=1, barrier=(2.77812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121282,'amu*angstrom^2'), symmetry=1, barrier=(2.78852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120963,'amu*angstrom^2'), symmetry=1, barrier=(2.78117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4214.4,'J/mol'), sigma=(6.97559,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=658.28 K, Pc=28.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28187,0.0729175,-0.000131969,1.35104e-07,-5.20697e-11,89235,30.5305], Tmin=(100,'K'), Tmax=(841.972,'K')), NASAPolynomial(coeffs=[1.35808,0.041739,-2.15235e-05,4.18441e-09,-2.89526e-13,90314.5,36.6625], Tmin=(841.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(741.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P) + radical(CCOJ) + radical(Cds_S) + radical(Cs_P)"""),
)

species(
    label = '[O][C]([O])CC1=CC1(20889)',
    structure = SMILES('[O][C]([O])CC1=CC1'),
    E0 = (509.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66945,0.0567312,-7.32772e-05,6.57253e-08,-2.54581e-11,61363.2,26.367], Tmin=(100,'K'), Tmax=(740.573,'K')), NASAPolynomial(coeffs=[4.16848,0.0363704,-1.81369e-05,3.57456e-09,-2.53321e-13,61181.3,16.3317], Tmin=(740.573,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(509.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(CCOJ) + radical(Cs_P) + radical(CCOJ)"""),
)

species(
    label = '[CH]=C1CC([O])([O])C1(20873)',
    structure = SMILES('[CH]=C1CC([O])([O])C1'),
    E0 = (426.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[12.8254,-0.0307266,0.000142121,-1.27755e-07,2.85057e-11,50870.5,-24.8688], Tmin=(100,'K'), Tmax=(1757.89,'K')), NASAPolynomial(coeffs=[78.4621,0.0190275,-7.02323e-05,1.72111e-08,-1.27424e-12,-2969.69,-466.057], Tmin=(1757.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Cds_P) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ)"""),
)

species(
    label = 'C=C1[CH]C([O])([O])C1(20890)',
    structure = SMILES('C=C1[CH]C([O])([O])C1'),
    E0 = (296.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[12.3536,-0.0266885,0.000139488,-1.2766e-07,2.87216e-11,35236.1,-25.6678], Tmin=(100,'K'), Tmax=(1748.01,'K')), NASAPolynomial(coeffs=[78.359,0.0203002,-7.07667e-05,1.73378e-08,-1.28512e-12,-18094.1,-467.509], Tmin=(1748.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]=[C]C[C]([O])[O](5460)',
    structure = SMILES('[CH]=[C]C[C]([O])[O]'),
    E0 = (765.014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,360,370,350,180,1978.08,1978.6],'cm^-1')),
        HinderedRotor(inertia=(0.146013,'amu*angstrom^2'), symmetry=1, barrier=(3.35713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146002,'amu*angstrom^2'), symmetry=1, barrier=(3.35688,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97061,0.0576059,-0.000116073,1.23952e-07,-4.83546e-11,92070.5,25.8144], Tmin=(100,'K'), Tmax=(852.847,'K')), NASAPolynomial(coeffs=[0.341154,0.0334913,-1.78048e-05,3.47479e-09,-2.39779e-13,93503.4,40.1871], Tmin=(852.847,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(765.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCOJ) + radical(CCOJ) + radical(Cds_S) + radical(Cs_P)"""),
)

species(
    label = '[C]C(=C)C[C]([O])[O](20891)',
    structure = SMILES('[C]C(=C)C[C]([O])[O]'),
    E0 = (910.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,360,370,350,180,180,1367.99,1368.1,1368.96],'cm^-1')),
        HinderedRotor(inertia=(0.276856,'amu*angstrom^2'), symmetry=1, barrier=(6.36547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.276549,'amu*angstrom^2'), symmetry=1, barrier=(6.3584,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23894,0.0709113,-0.000122872,1.19573e-07,-4.48332e-11,109597,26.774], Tmin=(100,'K'), Tmax=(828.418,'K')), NASAPolynomial(coeffs=[4.20752,0.0348124,-1.80984e-05,3.54181e-09,-2.4652e-13,109852,17.5188], Tmin=(828.418,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(910.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Cs_P) + radical(CCOJ) + radical(CCOJ) + radical(CJ3)"""),
)

species(
    label = 'N2',
    structure = SMILES('N#N'),
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'Ne',
    structure = SMILES('[Ne]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.1797,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ne""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'He',
    structure = SMILES('[He]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (4.0026,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(84.8076,'J/mol'), sigma=(2.576,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""He""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'Ar',
    structure = SMILES('[Ar]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (611.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (705.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (733.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (611.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (611.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (723.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (771.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (697.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (666.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1026.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (865.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1071.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (736.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (688.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (998.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (766.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (611.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (1028.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (843.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (764.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (814.281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (803.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (675.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (675.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (986.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (619.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (619.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (619.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (1180.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (1122.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]C(=C)C[C]([O])[O](18595)'],
    products = ['C=C([O])[O](1172)', 'C#C[CH2](17441)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]C(=C)C[C]([O])[O](18595)'],
    products = ['[CH]C1([CH2])CC1([O])[O](20876)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.17159e+11,'s^-1'), n=-0.0562706, Ea=(94.2611,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 91.8 to 94.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]C(=C)C[C]([O])[O](18595)'],
    products = ['[CH]C1([CH2])C[C]([O])O1(20824)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.448e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[CH]C(=C)[CH]C([O])=O(20877)'],
    products = ['[CH]C(=C)C[C]([O])[O](18595)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(84.475,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 80.3 to 84.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C([O])[O](1172)', '[CH][C]=C(18825)'],
    products = ['[CH]C(=C)C[C]([O])[O](18595)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.123469,'m^3/(mol*s)'), n=2.00579, Ea=(37.9107,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 36.0 to 37.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]C(=C)C[C]([O])[O](18595)'],
    products = ['[CH]C(=C)[CH]C([O])[O](20878)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.42705e+08,'s^-1'), n=1.50595, Ea=(111.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]C(=C)C[C]([O])[O](18595)'],
    products = ['[CH]C(=C)[CH][C]([O])O(20879)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.64e+07,'s^-1'), n=1.69, Ea=(159.41,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;O_rad_out;Cs_H_out_H/Cd] for rate rule [R3HJ;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]C(=[CH])CC([O])[O](20880)'],
    products = ['[CH]C(=C)C[C]([O])[O](18595)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]C(=[CH])C[C]([O])O(20881)'],
    products = ['[CH]C(=C)C[C]([O])[O](18595)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_3;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C]([O])[O](2304)', '[CH][C]=C(18825)'],
    products = ['[CH]C(=C)C[C]([O])[O](18595)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(8)', '[CH][C]([CH2])C=C([O])[O](20882)'],
    products = ['[CH]C(=C)C[C]([O])[O](18595)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(8)', '[CH]C(=[CH])C[C]([O])[O](20883)'],
    products = ['[CH]C(=C)C[C]([O])[O](18595)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.15742e+08,'m^3/(mol*s)'), n=0.0433333, Ea=(0.458029,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_rad;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]C(=C)C[C]([O])[O](18595)'],
    products = ['[CH][C]1CC([O])([O])C1(20884)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.78932e+07,'s^-1'), n=1.19089, Ea=(125.213,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_HH_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]C(=C)C[C]([O])[O](18595)'],
    products = ['[CH][C]1CO[C]([O])C1(20885)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.00945e+11,'s^-1'), n=0.246651, Ea=(76.6831,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra] for rate rule [R5_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]C([CH2])[CH][C]([O])[O](20886)'],
    products = ['[CH]C(=C)C[C]([O])[O](18595)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]C(=C)C([CH2])([O])[O](18614)'],
    products = ['[CH]C(=C)C[C]([O])[O](18595)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]C(=C)C[C]([O])[O](18595)'],
    products = ['[CH]C(=C)CC([O])=O(18773)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O][C][O](2319)', '[CH]C([CH2])=C(18776)'],
    products = ['[CH]C(=C)C[C]([O])[O](18595)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CH2(T)(28)', 'C#CC[C]([O])[O](5453)'],
    products = ['[CH]C(=C)C[C]([O])[O](18595)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(18.899,'m^3/(mol*s)'), n=1.76329, Ea=(16.1554,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-H;YJ] for rate rule [Ct-Cs_Ct-H;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][C]([O])[O](2304)', 'C#C[CH2](17441)'],
    products = ['[CH]C(=C)C[C]([O])[O](18595)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.0031273,'m^3/(mol*s)'), n=2.54618, Ea=(24.6367,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cdd_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]C(=C)C[C]([O])[O](18595)'],
    products = ['[CH][C](C)C=C([O])[O](20887)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.09894e+08,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]C(=C)C[C]([O])[O](18595)'],
    products = ['[CH2][C]([CH2])C=C([O])[O](18594)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]C(=C)C[C]([O])[O](18595)'],
    products = ['[CH]=C(C)[CH]C([O])=O(20888)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]C(=C)C[C]([O])[O](18595)'],
    products = ['[CH2]C(=C)[CH]C([O])=O(18587)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C]CC[C]([O])[O](18582)'],
    products = ['[CH]C(=C)C[C]([O])[O](18595)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]C(=C)C[C]([O])[O](18595)'],
    products = ['[O][C]([O])CC1=CC1(20889)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Ypri_rad_out] for rate rule [R3_SD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.82842712475
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]C(=C)C[C]([O])[O](18595)'],
    products = ['[CH]=C1CC([O])([O])C1(20873)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]C(=C)C[C]([O])[O](18595)'],
    products = ['C=C1[CH]C([O])([O])C1(20890)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(T)(28)', '[CH]=[C]C[C]([O])[O](5460)'],
    products = ['[CH]C(=C)C[C]([O])[O](18595)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(8)', '[C]C(=C)C[C]([O])[O](20891)'],
    products = ['[CH]C(=C)C[C]([O])[O](18595)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '4006',
    isomers = [
        '[CH]C(=C)C[C]([O])[O](18595)',
    ],
    reactants = [
        ('C=C([O])[O](1172)', 'C#C[CH2](17441)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4006',
    Tmin = (1200,'K'),
    Tmax = (1500,'K'),
    Tcount = 10,
    Tlist = ([1201.48,1213.22,1236.21,1269.31,1310.55,1356.92,1404.16,1447.02,1479.84,1497.7],'K'),
    Pmin = (1,'atm'),
    Pmax = (10,'atm'),
    Pcount = 10,
    Plist = ([1.02771,1.14872,1.41959,1.89986,2.67608,3.83649,5.40396,7.23219,8.93758,9.98989],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

