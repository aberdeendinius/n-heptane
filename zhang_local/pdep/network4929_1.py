species(
    label = '[CH]=[C]CC=C=C[O](22626)',
    structure = SMILES('[CH]=[C]CC=C=C[O]'),
    E0 = (646.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.06585,'amu*angstrom^2'), symmetry=1, barrier=(24.506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06786,'amu*angstrom^2'), symmetry=1, barrier=(24.5523,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.874198,0.0604722,-4.9312e-05,1.32472e-08,1.2772e-12,77928.4,27.0485], Tmin=(100,'K'), Tmax=(985.722,'K')), NASAPolynomial(coeffs=[16.5106,0.0141667,-4.93912e-06,8.83067e-10,-6.24589e-14,74012.8,-52.3884], Tmin=(985.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(646.936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C#CC=O(21959)',
    structure = SMILES('C#CC=O'),
    E0 = (84.2941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,346.083,346.148,1254.41,1569.49,1569.49],'cm^-1')),
        HinderedRotor(inertia=(0.289847,'amu*angstrom^2'), symmetry=1, barrier=(24.6472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3153.34,'J/mol'), sigma=(5.09602,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=492.54 K, Pc=54.07 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00307,0.020875,-1.61807e-05,6.30357e-09,-1.00066e-12,10174.9,9.76778], Tmin=(100,'K'), Tmax=(1462.04,'K')), NASAPolynomial(coeffs=[7.33456,0.00902438,-4.02236e-06,7.59526e-10,-5.26541e-14,8908.32,-12.7744], Tmin=(1462.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.2941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH)"""),
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
    label = '[CH]=C1CC1[C]=C[O](24784)',
    structure = SMILES('[CH]=C1CC1[C]=C[O]'),
    E0 = (676.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.981201,0.0563768,-3.4676e-05,-4.05219e-09,8.06236e-12,81452.4,22.7962], Tmin=(100,'K'), Tmax=(953.881,'K')), NASAPolynomial(coeffs=[17.0223,0.0129667,-3.92743e-06,6.76671e-10,-4.86827e-14,77306.8,-59.5244], Tmin=(953.881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(676.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_S) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = '[O]C=[C]C1C=[C]C1(24785)',
    structure = SMILES('[O]C=[C]C1C=[C]C1'),
    E0 = (639.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36876,0.0465842,-1.22713e-05,-2.28587e-08,1.3543e-11,76978.9,22.5998], Tmin=(100,'K'), Tmax=(965.255,'K')), NASAPolynomial(coeffs=[15.4347,0.0151074,-5.02228e-06,9.11612e-10,-6.66601e-14,73014.4,-51.2328], Tmin=(965.255,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(639.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Cds_S) + radical(C=COJ) + radical(cyclobutene-vinyl)"""),
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
    label = 'C#CC=C[C]=C[O](22635)',
    structure = SMILES('C#CC=C[C]=C[O]'),
    E0 = (454.145,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,750,770,3400,2100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2175,525,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.62614,'amu*angstrom^2'), symmetry=1, barrier=(37.3881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62728,'amu*angstrom^2'), symmetry=1, barrier=(37.4144,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.831519,0.0562972,-2.5861e-05,-2.38426e-08,1.8337e-11,54747.6,21.7645], Tmin=(100,'K'), Tmax=(911.449,'K')), NASAPolynomial(coeffs=[20.957,0.00392216,1.17357e-06,-3.44429e-10,2.26331e-14,49585.7,-81.6544], Tmin=(911.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = '[CH]=[C]CC=C=C=O(24786)',
    structure = SMILES('[CH]=[C]CC=C=C=O'),
    E0 = (691.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,1685,370,3010,987.5,1337.5,450,1655,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000525832,'amu*angstrom^2'), symmetry=1, barrier=(5.97033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000524487,'amu*angstrom^2'), symmetry=1, barrier=(5.95505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18904,0.0364105,-2.87925e-05,1.16452e-08,-1.89495e-12,83180,9.53559], Tmin=(100,'K'), Tmax=(1453.2,'K')), NASAPolynomial(coeffs=[10.3176,0.014036,-5.69729e-06,1.05001e-09,-7.22e-14,80817.5,-32.7183], Tmin=(1453.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(691.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=C[O](8556)',
    structure = SMILES('[CH]=C=C[O]'),
    E0 = (269.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.7401,0.0176211,1.62543e-05,-4.58526e-08,2.291e-11,32513.4,12.64], Tmin=(100,'K'), Tmax=(883.628,'K')), NASAPolynomial(coeffs=[13.5244,-0.00323064,4.17673e-06,-9.22668e-10,6.45049e-14,29515.7,-44.2318], Tmin=(883.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]C=CC=C=C[O](24787)',
    structure = SMILES('[CH]C=CC=C=C[O]'),
    E0 = (501.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.782723,0.0566767,-1.4346e-05,-3.18047e-08,1.91222e-11,60484.8,24.5867], Tmin=(100,'K'), Tmax=(940.73,'K')), NASAPolynomial(coeffs=[18.7528,0.0153334,-4.3366e-06,7.25523e-10,-5.27467e-14,55552.2,-69.2575], Tmin=(940.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet) + radical(C=COJ)"""),
)

species(
    label = '[CH]=[C]C[CH]C#CO(24788)',
    structure = SMILES('[CH]=[C]C[CH]C#CO'),
    E0 = (719.486,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,1685,370,2100,2250,500,550,328.605,329.542],'cm^-1')),
        HinderedRotor(inertia=(0.00158176,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107158,'amu*angstrom^2'), symmetry=1, barrier=(8.11403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.976918,'amu*angstrom^2'), symmetry=1, barrier=(74.5479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.976948,'amu*angstrom^2'), symmetry=1, barrier=(74.533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24856,0.062555,-7.31414e-05,4.58237e-08,-1.0733e-11,86631.5,26.823], Tmin=(100,'K'), Tmax=(762.294,'K')), NASAPolynomial(coeffs=[10.085,0.0230708,-8.99139e-06,1.5669e-09,-1.03497e-13,85084.3,-14.7212], Tmin=(762.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(719.486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CtH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(Cds_P) + radical(Sec_Propargyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC[C]=C=C[O](24789)',
    structure = SMILES('[CH]=CC[C]=C=C[O]'),
    E0 = (646.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.06757,'amu*angstrom^2'), symmetry=1, barrier=(24.5456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06396,'amu*angstrom^2'), symmetry=1, barrier=(24.4626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.874199,0.0604721,-4.93119e-05,1.32472e-08,1.2772e-12,77928.4,27.0485], Tmin=(100,'K'), Tmax=(985.722,'K')), NASAPolynomial(coeffs=[16.5106,0.0141667,-4.93912e-06,8.83067e-10,-6.2459e-14,74012.8,-52.3884], Tmin=(985.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(646.936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=C[C]=C[O](22614)',
    structure = SMILES('C=[C]C=C[C]=C[O]'),
    E0 = (476.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.63092,'amu*angstrom^2'), symmetry=1, barrier=(37.498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.628,'amu*angstrom^2'), symmetry=1, barrier=(37.4309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.369045,0.0747065,-7.9467e-05,4.04014e-08,-7.54556e-12,57467.9,26.7793], Tmin=(100,'K'), Tmax=(1559.56,'K')), NASAPolynomial(coeffs=[20.0548,0.00664978,1.06533e-06,-4.6779e-10,3.88295e-14,53003.4,-74.7195], Tmin=(1559.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C[C]=C=C[O](24306)',
    structure = SMILES('C=[C]C[C]=C=C[O]'),
    E0 = (637.682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.996479,'amu*angstrom^2'), symmetry=1, barrier=(22.911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.99986,'amu*angstrom^2'), symmetry=1, barrier=(22.9888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.788263,0.0629847,-6.2157e-05,3.07557e-08,-5.93811e-12,76817.5,27.4308], Tmin=(100,'K'), Tmax=(1268.94,'K')), NASAPolynomial(coeffs=[16.0234,0.0149606,-5.38889e-06,9.31594e-10,-6.23959e-14,72951,-49.6988], Tmin=(1268.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(637.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]CC#C[CH]O(24790)',
    structure = SMILES('[CH]=[C]CC#C[CH]O'),
    E0 = (728.049,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,1685,370,2100,2250,500,550,301.461,1136.91],'cm^-1')),
        HinderedRotor(inertia=(0.147023,'amu*angstrom^2'), symmetry=1, barrier=(9.48199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147022,'amu*angstrom^2'), symmetry=1, barrier=(9.48186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147027,'amu*angstrom^2'), symmetry=1, barrier=(9.48198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00323,'amu*angstrom^2'), symmetry=1, barrier=(64.6982,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00611,0.0697097,-9.97927e-05,7.85423e-08,-2.4216e-11,87668.2,28.3382], Tmin=(100,'K'), Tmax=(904.663,'K')), NASAPolynomial(coeffs=[9.50469,0.0236755,-9.44136e-06,1.62638e-09,-1.04813e-13,86476.6,-9.89822], Tmin=(904.663,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(728.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CtHH) + group(Cs-CtOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Cds_P) + radical(Cds_S) + radical(CCsJOH)"""),
)

species(
    label = '[CH]=CC[CH][C]=C=O(24791)',
    structure = SMILES('[CH]=CC[CH][C]=C=O'),
    E0 = (575.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43625,0.0594023,-6.06644e-05,3.43828e-08,-8.14418e-12,69272.5,25.8384], Tmin=(100,'K'), Tmax=(1001.34,'K')), NASAPolynomial(coeffs=[9.44954,0.0273921,-1.27136e-05,2.45838e-09,-1.73795e-13,67667.7,-12.8318], Tmin=(1001.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(575.217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(CCJC(C)=C=O) + radical(CCCJ=C=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=C[CH][C]=CO(24792)',
    structure = SMILES('[CH]=C=C[CH][C]=CO'),
    E0 = (514.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.892675,0.0521838,-4.83939e-06,-5.17907e-08,3.04716e-11,62016.8,25.4694], Tmin=(100,'K'), Tmax=(889.314,'K')), NASAPolynomial(coeffs=[22.1894,0.00172589,3.80739e-06,-9.54788e-10,6.78885e-14,56436.3,-84.8552], Tmin=(889.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(514.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=CCJC=C) + radical(C=C=CJ)"""),
)

species(
    label = 'C=[C]C[CH][C]=C=O(24307)',
    structure = SMILES('C=[C]C[CH][C]=C=O'),
    E0 = (565.963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2120,512.5,787.5,3025,407.5,1350,352.5,304.114,305.404,1203.3],'cm^-1')),
        HinderedRotor(inertia=(0.159997,'amu*angstrom^2'), symmetry=1, barrier=(10.5335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.886339,'amu*angstrom^2'), symmetry=1, barrier=(57.9697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307049,'amu*angstrom^2'), symmetry=1, barrier=(20.1994,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42967,0.0609813,-7.03206e-05,4.80318e-08,-1.39051e-11,68158.3,25.9364], Tmin=(100,'K'), Tmax=(825.376,'K')), NASAPolynomial(coeffs=[7.79175,0.0301483,-1.42851e-05,2.7704e-09,-1.95504e-13,67108.1,-3.53579], Tmin=(825.376,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCJC(C)=C=O) + radical(CCCJ=C=O)"""),
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
    label = '[CH]=C=C[CH][C]=C[O](22638)',
    structure = SMILES('[CH]=C=C[CH][C]=C[O]'),
    E0 = (656.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.84007,'amu*angstrom^2'), symmetry=1, barrier=(42.3067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.84024,'amu*angstrom^2'), symmetry=1, barrier=(42.3107,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30798,0.0451968,-1.54147e-06,-4.51244e-08,2.53584e-11,79013.7,25.3872], Tmin=(100,'K'), Tmax=(902.276,'K')), NASAPolynomial(coeffs=[18.8022,0.00563787,1.05541e-06,-3.69868e-10,2.58024e-14,74310.1,-65.784], Tmin=(902.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(656.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C=CJ) + radical(C=CCJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C[C]=C=C[O](24793)',
    structure = SMILES('[CH]=[C]C[C]=C=C[O]'),
    E0 = (884.778,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.02673,'amu*angstrom^2'), symmetry=1, barrier=(23.6066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02283,'amu*angstrom^2'), symmetry=1, barrier=(23.5169,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.874011,0.0646744,-7.26473e-05,4.05177e-08,-8.77503e-12,106530,27.5235], Tmin=(100,'K'), Tmax=(1136.36,'K')), NASAPolynomial(coeffs=[15.4021,0.0135351,-5.1432e-06,9.15066e-10,-6.23981e-14,103228,-44.4233], Tmin=(1136.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(884.778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C[CH][C]=C=O(24794)',
    structure = SMILES('[CH]=[C]C[CH][C]=C=O'),
    E0 = (813.059,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,3120,650,792.5,1650,3025,407.5,1350,352.5,292.55,293.789],'cm^-1')),
        HinderedRotor(inertia=(0.0688635,'amu*angstrom^2'), symmetry=1, barrier=(4.11457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.355583,'amu*angstrom^2'), symmetry=1, barrier=(21.6174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.914852,'amu*angstrom^2'), symmetry=1, barrier=(55.5995,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18871,0.0666615,-9.54696e-05,7.79543e-08,-2.59466e-11,97885.1,27.1924], Tmin=(100,'K'), Tmax=(771.998,'K')), NASAPolynomial(coeffs=[8.45671,0.0265588,-1.28001e-05,2.4625e-09,-1.71421e-13,96835.8,-5.51889], Tmin=(771.998,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(813.059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(CCCJ=C=O) + radical(Cds_S) + radical(CCJC(C)=C=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC=C1[CH]O1(24795)',
    structure = SMILES('[CH]=[C]CC=C1[CH]O1'),
    E0 = (673.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28803,0.0393593,2.76673e-05,-7.8929e-08,3.71186e-11,81073.8,24.1447], Tmin=(100,'K'), Tmax=(932.221,'K')), NASAPolynomial(coeffs=[21.6166,0.00400391,1.09274e-06,-2.36694e-10,8.25534e-15,75029.8,-84.5907], Tmin=(932.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(673.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(Cds_S) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C[CH]C1=C[O](24796)',
    structure = SMILES('[CH]=C1C[CH]C1=C[O]'),
    E0 = (507.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5572,0.0335272,4.13133e-05,-8.82117e-08,3.90565e-11,61089,19.9963], Tmin=(100,'K'), Tmax=(938.183,'K')), NASAPolynomial(coeffs=[19.3873,0.00878142,-1.10115e-06,1.81323e-10,-2.0857e-14,55486.9,-76.9121], Tmin=(938.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(507.037,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(Cds_P) + radical(C=COJ) + radical(Allyl_S)"""),
)

species(
    label = '[CH]=[C]CC1[C]=CO1(24797)',
    structure = SMILES('[CH]=[C]CC1[C]=CO1'),
    E0 = (765.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06504,0.0448609,1.52652e-05,-6.78365e-08,3.34898e-11,92192.8,24.6248], Tmin=(100,'K'), Tmax=(934.158,'K')), NASAPolynomial(coeffs=[22.444,0.00349791,1.10726e-06,-2.29581e-10,7.73759e-15,86009,-88.7792], Tmin=(934.158,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(765.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[O]C=C1[CH]C[C]=C1(24798)',
    structure = SMILES('[O]C=C1[CH]C[C]=C1'),
    E0 = (428.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87203,0.0252799,6.10855e-05,-1.05735e-07,4.448e-11,51602.4,19.5098], Tmin=(100,'K'), Tmax=(941.537,'K')), NASAPolynomial(coeffs=[18.2956,0.0101021,-1.71259e-06,3.15951e-10,-3.16067e-14,46089.8,-71.5861], Tmin=(941.537,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(428.246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(3-Methylenecyclopentene) + radical(cyclopentene-vinyl) + radical(Allyl_S) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C1CC=[C]C1[O](24799)',
    structure = SMILES('[CH]=C1CC=[C]C1[O]'),
    E0 = (672.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95376,0.0336542,1.21327e-05,-3.69792e-08,1.53496e-11,80966.7,24.4191], Tmin=(100,'K'), Tmax=(1029.21,'K')), NASAPolynomial(coeffs=[11.7866,0.021871,-9.21664e-06,1.80252e-09,-1.31803e-13,77542.8,-30.1028], Tmin=(1029.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(672.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(4-Methylenecyclopentene) + radical(cyclopentene-vinyl) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[O]C1[C]=CC[C]=C1(24800)',
    structure = SMILES('[O]C1[C]=CC[C]=C1'),
    E0 = (622.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8934,0.0353991,7.96737e-06,-3.30022e-08,1.40555e-11,74898.1,23.3562], Tmin=(100,'K'), Tmax=(1030.23,'K')), NASAPolynomial(coeffs=[11.7904,0.0220907,-9.22725e-06,1.79025e-09,-1.30179e-13,71525.9,-31.1551], Tmin=(1030.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(622.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cyclohexadiene) + radical(Cds_S) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C#CC=C[C]=CO(24801)',
    structure = SMILES('C#CC=C[C]=CO'),
    E0 = (312.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.01417,0.0798868,-8.62717e-05,4.27272e-08,-7.6394e-12,37813.8,26.9993], Tmin=(100,'K'), Tmax=(1664.17,'K')), NASAPolynomial(coeffs=[22.6428,0.00174304,3.34568e-06,-8.5801e-10,6.25937e-14,32886.9,-90.3268], Tmin=(1664.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=CCC=C=C=O(24802)',
    structure = SMILES('[CH]=CCC=C=C=O'),
    E0 = (453.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19091,0.0322241,-5.71624e-06,-1.49453e-08,7.70619e-12,54578,9.05264], Tmin=(100,'K'), Tmax=(1024.77,'K')), NASAPolynomial(coeffs=[10.6207,0.0159305,-6.17988e-06,1.1733e-09,-8.47214e-14,51978.1,-36.0779], Tmin=(1024.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P)"""),
)

species(
    label = '[CH][C]=CC[C]=C[O](21989)',
    structure = SMILES('[CH][C]=CC[C]=C[O]'),
    E0 = (831.762,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,262.359,262.567,263.012,263.317,263.547,263.804],'cm^-1')),
        HinderedRotor(inertia=(1.00729,'amu*angstrom^2'), symmetry=1, barrier=(49.6533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00784,'amu*angstrom^2'), symmetry=1, barrier=(49.6507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00962,'amu*angstrom^2'), symmetry=1, barrier=(49.6696,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.675442,0.0640658,-5.52067e-05,2.4507e-08,-4.33321e-12,100165,29.8951], Tmin=(100,'K'), Tmax=(1364.77,'K')), NASAPolynomial(coeffs=[15.4049,0.0208953,-7.75853e-06,1.32937e-09,-8.75057e-14,96144.7,-45.7466], Tmin=(1364.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(831.762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=COJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C[C]C=C[O](24803)',
    structure = SMILES('[CH]=[C]C[C]C=C[O]'),
    E0 = (929.523,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,201.821,201.912,201.915,201.926],'cm^-1')),
        HinderedRotor(inertia=(0.781012,'amu*angstrom^2'), symmetry=1, barrier=(22.5956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.780664,'amu*angstrom^2'), symmetry=1, barrier=(22.5957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.780808,'amu*angstrom^2'), symmetry=1, barrier=(22.5957,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.723389,0.0639226,-5.62094e-05,1.85316e-08,-1.70775e-13,111921,28.8529], Tmin=(100,'K'), Tmax=(983.283,'K')), NASAPolynomial(coeffs=[17.2457,0.0134055,-4.61495e-06,8.18883e-10,-5.77822e-14,107865,-54.6835], Tmin=(983.283,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(929.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(C=COJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C[CH]C=[C][O](24804)',
    structure = SMILES('[CH]=[C]C[CH]C=[C][O]'),
    E0 = (862.165,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,3025,407.5,1350,352.5,866.163,868.918,877.583],'cm^-1')),
        HinderedRotor(inertia=(0.524294,'amu*angstrom^2'), symmetry=1, barrier=(12.0546,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.52316,'amu*angstrom^2'), symmetry=1, barrier=(12.0285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208303,'amu*angstrom^2'), symmetry=1, barrier=(27.8512,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26559,0.0563946,-5.08397e-05,2.37311e-08,-4.45327e-12,103796,29.0254], Tmin=(100,'K'), Tmax=(1276.54,'K')), NASAPolynomial(coeffs=[12.812,0.0202145,-8.32645e-06,1.52892e-09,-1.05171e-13,100848,-29.4985], Tmin=(1276.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(862.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Allyl_S) + radical(Cds_P) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = '[CH][C]=C[CH]C=C[O](21988)',
    structure = SMILES('[CH][C]=C[CH]C=C[O]'),
    E0 = (695.645,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,343.609,343.867,343.935,344.076,344.096,344.323],'cm^-1')),
        HinderedRotor(inertia=(0.593938,'amu*angstrom^2'), symmetry=1, barrier=(49.887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.595243,'amu*angstrom^2'), symmetry=1, barrier=(49.8834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.5945,'amu*angstrom^2'), symmetry=1, barrier=(49.8854,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25666,0.0444989,1.5576e-05,-5.93053e-08,2.81725e-11,83780.1,26.9983], Tmin=(100,'K'), Tmax=(936.138,'K')), NASAPolynomial(coeffs=[16.9984,0.0172117,-4.75425e-06,7.88015e-10,-5.73994e-14,79081.2,-57.2631], Tmin=(936.138,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(695.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(C=CCJC=C) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C=C[CH][C]C[O](24805)',
    structure = SMILES('[CH]=C=C[CH][C]C[O]'),
    E0 = (909.959,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,180,180,180,1153.69],'cm^-1')),
        HinderedRotor(inertia=(0.00803636,'amu*angstrom^2'), symmetry=1, barrier=(62.9046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.73379,'amu*angstrom^2'), symmetry=1, barrier=(62.8551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.73434,'amu*angstrom^2'), symmetry=1, barrier=(62.8679,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20796,0.0659135,-8.81524e-05,7.10541e-08,-2.35617e-11,109539,26.1952], Tmin=(100,'K'), Tmax=(809.486,'K')), NASAPolynomial(coeffs=[7.40646,0.0305372,-1.3803e-05,2.57798e-09,-1.76421e-13,108691,-1.43832], Tmin=(809.486,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(909.959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_S) + radical(C=C=CJ) + radical(CCOJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)C=C=C[O](22618)',
    structure = SMILES('[CH]C(=C)C=C=C[O]'),
    E0 = (500.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.00191,'amu*angstrom^2'), symmetry=1, barrier=(46.0278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00636,'amu*angstrom^2'), symmetry=1, barrier=(46.1302,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.61675,0.0607309,-2.30026e-05,-2.49436e-08,1.72533e-11,60302.2,23.9773], Tmin=(100,'K'), Tmax=(936.072,'K')), NASAPolynomial(coeffs=[19.4917,0.0144273,-3.85094e-06,6.21069e-10,-4.4817e-14,55263.5,-73.8753], Tmin=(936.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(500.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]=[CH](18830)',
    structure = SMILES('[C]=[CH]'),
    E0 = (847.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([562.459,1392.74,3112.83],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86574,0.0023643,1.31489e-06,-2.64796e-09,1.00932e-12,101918,6.62295], Tmin=(100,'K'), Tmax=(1016.96,'K')), NASAPolynomial(coeffs=[4.17915,0.00251347,-9.43473e-07,1.6872e-10,-1.15831e-14,101782,4.75427], Tmin=(1016.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'C=C[C]=C[O](22438)',
    structure = SMILES('C=C[C]=C[O]'),
    E0 = (225.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.61747,'amu*angstrom^2'), symmetry=1, barrier=(37.1889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98077,0.0322088,7.53766e-06,-4.54875e-08,2.38742e-11,27216,16.1107], Tmin=(100,'K'), Tmax=(899.949,'K')), NASAPolynomial(coeffs=[16.1069,0.00249504,1.93925e-06,-5.05324e-10,3.47549e-14,23334.1,-57.9915], Tmin=(899.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = 'O(T)(63)',
    structure = SMILES('[O]'),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,29230.2,4.09104], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,29230.2,4.09104], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=[C]C[CH]C#C(19875)',
    structure = SMILES('[CH]=[C]C[CH]C#C'),
    E0 = (860.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,750,770,3400,2100,3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,1685,370,2175,525,299.094],'cm^-1')),
        HinderedRotor(inertia=(0.105027,'amu*angstrom^2'), symmetry=1, barrier=(6.60885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1718,'amu*angstrom^2'), symmetry=1, barrier=(74.1119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16982,'amu*angstrom^2'), symmetry=1, barrier=(74.1037,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63213,0.0531178,-5.49663e-05,2.62977e-08,-2.3909e-12,103633,22.1139], Tmin=(100,'K'), Tmax=(742.746,'K')), NASAPolynomial(coeffs=[9.87547,0.0188677,-6.28284e-06,9.88323e-10,-6.10194e-14,102129,-17.0874], Tmin=(742.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(860.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(Cds_P) + radical(Sec_Propargyl)"""),
)

species(
    label = '[C]=[C]CC=C=C[O](24806)',
    structure = SMILES('[C]=[C]CC=C=C[O]'),
    E0 = (957.942,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.994823,'amu*angstrom^2'), symmetry=1, barrier=(22.8729,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.999535,'amu*angstrom^2'), symmetry=1, barrier=(22.9813,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.954221,0.0627962,-6.78491e-05,3.64361e-08,-7.63102e-12,115327,26.7374], Tmin=(100,'K'), Tmax=(1170.65,'K')), NASAPolynomial(coeffs=[15.1503,0.0142891,-5.69467e-06,1.03992e-09,-7.18934e-14,112003,-43.9872], Tmin=(1170.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(957.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=[C]CC#CC=O(24807)',
    structure = SMILES('[CH]=[C]CC#CC=O'),
    E0 = (613.107,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2100,2250,500,550,3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,260.199],'cm^-1')),
        HinderedRotor(inertia=(0.278896,'amu*angstrom^2'), symmetry=1, barrier=(13.4002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.278915,'amu*angstrom^2'), symmetry=1, barrier=(13.4001,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.747405,'amu*angstrom^2'), symmetry=1, barrier=(35.9072,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74121,0.0536148,-5.92387e-05,3.81123e-08,-1.04039e-11,73817.7,23.8158], Tmin=(100,'K'), Tmax=(871.145,'K')), NASAPolynomial(coeffs=[7.63181,0.0265668,-1.26648e-05,2.46987e-09,-1.75097e-13,72791.4,-3.79027], Tmin=(871.145,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(613.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CtHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C[C]=CC=O(24808)',
    structure = SMILES('[CH]=[C]C[C]=CC=O'),
    E0 = (688.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.529181,'amu*angstrom^2'), symmetry=1, barrier=(12.1669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.5288,'amu*angstrom^2'), symmetry=1, barrier=(12.1581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.528691,'amu*angstrom^2'), symmetry=1, barrier=(12.1556,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29346,0.0655943,-9.31412e-05,8.16507e-08,-3.01388e-11,82844,26.9657], Tmin=(100,'K'), Tmax=(722.396,'K')), NASAPolynomial(coeffs=[6.35106,0.0333218,-1.72679e-05,3.45232e-09,-2.46352e-13,82224.6,4.98115], Tmin=(722.396,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(688.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C[CH]C=C=O(24809)',
    structure = SMILES('[CH]=[C]C[CH]C=C=O'),
    E0 = (610.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44097,0.0622572,-7.31495e-05,4.49817e-08,-8.59548e-12,73546.3,26.61], Tmin=(100,'K'), Tmax=(615.069,'K')), NASAPolynomial(coeffs=[7.65437,0.0299968,-1.43436e-05,2.77893e-09,-1.95392e-13,72627.8,-1.59901], Tmin=(615.069,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(610.778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(CCJC(C)=C=O) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C=C[CH]C=C[O](21998)',
    structure = SMILES('[CH]=C=C[CH]C=C[O]'),
    E0 = (418.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.81308,'amu*angstrom^2'), symmetry=1, barrier=(41.6863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.81229,'amu*angstrom^2'), symmetry=1, barrier=(41.668,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33213,0.0407037,2.28675e-05,-7.3903e-08,3.61183e-11,50410.7,24.8266], Tmin=(100,'K'), Tmax=(908.639,'K')), NASAPolynomial(coeffs=[20.0131,0.00610111,1.35435e-06,-4.23921e-10,2.75502e-14,45049.5,-74.3289], Tmin=(908.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CCJC=C) + radical(C=COJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C1CC=[C][CH]O1(24810)',
    structure = SMILES('[CH]=C1CC=[C][CH]O1'),
    E0 = (519.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14585,0.0209352,6.63262e-05,-9.85149e-08,3.75416e-11,62557.2,18.7142], Tmin=(100,'K'), Tmax=(997.291,'K')), NASAPolynomial(coeffs=[14.3743,0.0215348,-9.2473e-06,1.92061e-09,-1.48444e-13,57649.2,-52.6258], Tmin=(997.291,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(519.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(C=CCJ(O)C) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[C]1[CH]OC=[C]CC=1(24811)',
    structure = SMILES('[C]1[CH]OC=[C]CC=1'),
    E0 = (546.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04399,-0.0153176,0.000246569,-3.57549e-07,1.50985e-10,65807.2,22.6827], Tmin=(100,'K'), Tmax=(899.037,'K')), NASAPolynomial(coeffs=[41.4069,-0.0372165,2.74413e-05,-5.47461e-09,3.63508e-13,52536.7,-197.472], Tmin=(899.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(Cds_S) + radical(Cds_S) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C#CC=CC=C[O](24812)',
    structure = SMILES('C#CC=CC=C[O]'),
    E0 = (255.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.904079,0.0498327,2.8764e-06,-5.61197e-08,2.97761e-11,30815.7,21.9065], Tmin=(100,'K'), Tmax=(927.719,'K')), NASAPolynomial(coeffs=[22.4567,0.00390139,1.15485e-06,-2.78096e-10,1.33117e-14,24794.4,-91.3554], Tmin=(927.719,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C1CC=C1C=O(24813)',
    structure = SMILES('[CH]=C1CC=C1C=O'),
    E0 = (352.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.1033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29701,0.0420962,1.28977e-05,-5.55565e-08,2.61108e-11,42553.1,19.8563], Tmin=(100,'K'), Tmax=(968.145,'K')), NASAPolynomial(coeffs=[19.2194,0.0103077,-3.32641e-06,7.02236e-10,-5.91857e-14,37102.2,-76.2573], Tmin=(968.145,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene) + radical(Cds_P)"""),
)

species(
    label = '[CH]=O(373)',
    structure = SMILES('[CH]=O'),
    E0 = (33.3613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([880.033,2045.55,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.35602,-0.00347087,1.25664e-05,-9.99488e-09,2.27888e-12,3995.77,2.75113], Tmin=(100,'K'), Tmax=(1565.72,'K')), NASAPolynomial(coeffs=[4.61867,0.00504456,-4.39241e-06,9.73282e-10,-7.07436e-14,2787.5,-2.22961], Tmin=(1565.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.3613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-OdHH) + radical(HCdsJO)"""),
)

species(
    label = '[C]=CC[C]=[CH](19256)',
    structure = SMILES('[C]=CC[C]=[CH]'),
    E0 = (1131.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,238.314],'cm^-1')),
        HinderedRotor(inertia=(0.266204,'amu*angstrom^2'), symmetry=1, barrier=(10.7291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266216,'amu*angstrom^2'), symmetry=1, barrier=(10.7291,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30756,0.0399091,-4.4945e-05,2.98958e-08,-8.40632e-12,136181,20.1622], Tmin=(100,'K'), Tmax=(849.582,'K')), NASAPolynomial(coeffs=[6.6232,0.0195903,-9.07073e-06,1.74536e-09,-1.22727e-13,135448,0.0451893], Tmin=(849.582,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1131.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(CdCdJ2_triplet)"""),
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
    E0 = (646.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (730.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (700.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (677.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (905.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (646.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (764.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (897.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (788.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (691.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1042.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1043.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (679.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (737.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (798.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (884.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (867.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (1096.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (1024.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (834.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (773.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (772.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (679.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (716.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (659.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (671.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (671.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (859.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (951.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (885.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (759.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (934.932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (892.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (1107.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (1267.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (1169.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (834.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (715.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (884.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (809.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (823.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (705.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (737.733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (735.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (655.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (1199.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['C#CC=O(21959)', 'C#C[CH2](17441)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[CH]=C1CC1[C]=C[O](24784)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.8958e+11,'s^-1'), n=-0.055489, Ea=(83.6851,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[O]C=[C]C1C=[C]C1(24785)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C#CC=C[C]=C[O](22635)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(197.535,'m^3/(mol*s)'), n=1.6075, Ea=(11.7084,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cdd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', '[CH]=[C]CC=C=C=O(24786)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C=C[O](8556)', 'C#C[CH2](17441)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(8.86332e-05,'m^3/(mol*s)'), n=2.74069, Ea=(48.5809,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cd_R;CdsJ=Cdd] + [Cds_Cdd;CJ] for rate rule [Cds_Cdd;CdsJ=Cdd]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 43.6 to 48.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[CH]C=CC=C=C[O](24787)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.18474e+10,'s^-1'), n=0.962222, Ea=(117.477,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]C[CH]C#CO(24788)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_double;XH_out] for rate rule [R2H_S;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=CC[C]=C=C[O](24789)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_double;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['C=[C]C=C[C]=C[O](22614)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/Cd] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['C=[C]C[C]=C=C[O](24306)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.15e+11,'s^-1'), n=0.388, Ea=(395.848,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R4HJ_1;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]CC#C[CH]O(24790)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;XH_out] for rate rule [R4H_MMS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[CH]=CC[CH][C]=C=O(24791)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[CH]=C=C[CH][C]=CO(24792)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(492144,'s^-1'), n=1.73305, Ea=(91.0472,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;O_rad_out;XH_out] for rate rule [R5H_SMMS;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['C=[C]C[CH][C]=C=O(24307)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C=C[O](8556)', '[CH][C]=C(18825)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(40577.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Cd_allenic]
Euclidian distance = 0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', '[CH]=C=C[CH][C]=C[O](22638)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[CH]=[C]C[C]=C=C[O](24793)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[CH]=[C]C[CH][C]=C=O(24794)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[CH]=[C]CC=C1[CH]O1(24795)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[CH]=C1C[CH]C1=C[O](24796)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[CH]=[C]CC1[C]=CO1(24797)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.47276e+11,'s^-1'), n=0.183155, Ea=(125.658,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[O]C=C1[CH]C[C]=C1(24798)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(9.13977e+10,'s^-1'), n=0.0396934, Ea=(32.3534,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[CH]=C1CC=[C]C1[O](24799)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.01304e+12,'s^-1'), n=-0.3725, Ea=(69.6427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[O]C1[C]=CC[C]=C1(24800)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.16959e+10,'s^-1'), n=0.31, Ea=(12.393,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra;radadd_intra_cdsingleH] for rate rule [R6_linear;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['C#CC=C[C]=CO(24801)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[CH]=CCC=C=C=O(24802)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH][C]=CC[C]=C[O](21989)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=[C]C[C]C=C[O](24803)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=[C]C[CH]C=[C][O](24804)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH][C]=C[CH]C=C[O](21988)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C=C[CH][C]C[O](24805)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[CH]C(=C)C=C=C[O](22618)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[C]=[CH](18830)', 'C=C[C]=C[O](22438)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O(T)(63)', '[CH]=[C]C[CH]C#C(19875)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(8)', '[C]=[C]CC=C=C[O](24806)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(8)', '[CH]=[C]CC#CC=O(24807)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3249.46,'m^3/(mol*s)'), n=1.38433, Ea=(9.80868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-De;HJ] for rate rule [Ct-Cs_Ct-CO;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C#CC=O(21959)', '[CH][C]=C(18825)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(0.106247,'m^3/(mol*s)'), n=2.32278, Ea=(16.475,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-De;CJ] for rate rule [Ct-H_Ct-CO;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=[C]C[C]=CC=O(24808)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[CH]=[C]C[CH]C=C=O(24809)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[CH]=C=C[CH]C=C[O](21998)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.65652e+07,'s^-1'), n=1.67955, Ea=(176.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[CH]=C1CC=[C][CH]O1(24810)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(4.39512e+11,'s^-1'), n=0.277081, Ea=(58.6284,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[C]1[CH]OC=[C]CC=1(24811)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(3.41957e+11,'s^-1'), n=0.2055, Ea=(90.797,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7_linear;multiplebond_intra;radadd_intra_cdsingleH] + [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['C#CC=CC=C[O](24812)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=[C]CC=C=C[O](22626)'],
    products = ['[CH]=C1CC=C1C=O(24813)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=O(373)', '[C]=CC[C]=[CH](19256)'],
    products = ['[CH]=[C]CC=C=C[O](22626)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '4929',
    isomers = [
        '[CH]=[C]CC=C=C[O](22626)',
    ],
    reactants = [
        ('C#CC=O(21959)', 'C#C[CH2](17441)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4929',
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

