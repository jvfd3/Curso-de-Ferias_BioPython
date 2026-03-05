
/*
//
// peptide.scad
// Copyright (c) 2019 Robert T. Miller.  All rights reserved.
// This file is part of the Biopython distribution and governed by your
// choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
// Please see the LICENSE file that should have been included as part of this
// package.
//
// This is the support file to build an OpenSCAD (http://www.openscad.org/) model
// of a protein from internal coordinates.  The resulting model may be constructed
// on a 3D printer.
//
// data matrices should be appended below to form a program ready to load into
// the OpenSCAD application.
//
//  The protein_scale value used throughout is the second element of the
//    protein[] array appended below.
//    This is the value supplied when generating the data for build units per
//    PDB angstrom.
//    You may wish to modify it here to adjust the appearance of the model in
//    terms of atom sphere or bond cylinder diameter, however the bond lengths
//    are fixed with the supplied value when the data matrices are generated.
//    Atom sphere and bond cylinder radii may be individually adjusted below as
//    well.
//
//  $fn (fragment number) is an OpenSCAD parameter controlling the smoothness
//    of the model surface.  Smaller values will render faster, but yield more
//    'blocky' models.
//
//  This is intended to be a working example, you are encouraged to modify the
//    OpenSCAD subroutines below to generate a model to your liking.  For more
//    information, start with http://www.openscad.org/cheatsheet/index.html
//
//  Note especially the hedronDispatch() subroutine below: here you may select
//    hedra based on residue, sequence position, and class (hedron atoms) for
//    special handling.  Also see the per hedron render options in the hedra[]
//    array.
//
//  Use the colorSel variable to restrict rendered objects by assigned color,
       and thereby separate STL files for a multi-material printer.

//  If you modify this file, you may find it useful to generate the data
//    matrices without this OpenSCAD code by calling write_SCAD() with the
//    includeCode=False option, then use the OpenSCAD 'include <>' facility at
//    the end of your modified OpenSCAD program.
*/

rotate([-90,0,0])  // convenient for default location (no N-Ca-C start coordinates)
    chain(protein);   // this is the main subroutine call to  build the structure

// top-level OpenSCAD $fn for visible surfaces.  Rotatable bonds use $fn=8
// inside, regardless of this setting.

$fn = 20;  // $fn=8 should print with minimal support

tubes=false;     // style: render atoms and bonds as constant diameter cylinders, preferred for rotatable bonds / h-bonds
support=false;   // enable print-in-place internal support for rotatable bonds
// N.B. rotatable bonds must be parallel to build plate for internal support
// structures to be generated correctly by slicer

// output parameters:

atomScale=1.0;  // 0.8 better for rotatable bonds
defaultAtomRadius = 0.77;  // used if tubes = true

bondRadius = (tubes ? defaultAtomRadius * atomScale : 0.4);
jBondRadius = defaultAtomRadius * atomScale;  // radius for rotatable bonds

// general printer, slicer, print settings:

layerHeight=0.15;  // must match slicer setting for print-in-place support
clearance=0.3;     // sliding clearance - can be smaller (0.2) if not doing print-in-place
pClearance=0.2;    // press-fit clearance (magnets for h-bonds)
shim=0.05;         // extra to make OpenSCAD surfaces distinct in difference()
nozzleDiameter=0.4;

// need one magnet for each side of hydrogen bond, suggest 3mm x 5mm e.g. from eBay
// use compass to identify poles if you care, North pointing (red) repels compass North pointer
magR=3/2;    // magnet radius
magL=5;      // magnet length

// for $fn=8 which works nice on fdm printer:

oRot = 22.5;              // 45/2, rotate to make fn=8 spheres and cylinders flat on build plate
apmFac = cos(180/8);      // apothem factor - multiply by radius for center to octagon side distance
octSide = 2* tan(180/8);  // multiply by radius to get length of octagon side
// for values of $fn:
fnRot = ($fn ? 90-(180/$fn) : 90-(180/30));

bondLenFac = 0.6;         // fraction of bond length to extend from atom for each arm of hedron in join

hblen = 1.97;             // hydrogen bond length

wall = 3*nozzleDiameter;
joinerStep = 1;           // radius difference between rotatable bond axle and end knob inside bond cylinder

caTop = false;     // only make top of N_C-alpha_C hedron plus C-beta (see hedron() and hedron_dispatch() examples)

colorSel = "all";  // [all, red, green, blue, gray, yellow, black, white]
bondColor = "white";  // [black, white, green]

/*
//
// Selector to restrict output according to assigned color.
// See https://nextjeff.com/creating-multi-extruder-designs-in-openscad-for-3d-printing-6c43a002ef64
//
*/
module colorSelect(Color) {
    color(Color) {
        if (colorSel == "all" || colorSel == Color) {
            children();
        }
    }
}

/*
//
// Generate a sphere to represent an atom.
// Colour and size determined for the atom covalent radius specified by the
//   parameter 'a' by lookup in the atomData table below, then scaled by the
//   supplied parameter 'scal'.
//
// scal : protein_scale
// clr : additional radius if used to create clearance for rotatable bonds
//
*/
module atom(a,scal,clr=0)
{
    ad = atomData[search([a],atomData)[0]];
    colorSelect(ad[1]) {
        rotate([0,0,fnRot]) sphere(r=((ad[2]*atomScale)*scal)+clr);
    }
}

/*
//
// a hedron (below) may be 'reversed' in terms of the order of its two bonds;
// this function fixes the ordering
//
*/
function hFlip(h,rev) =
        //   yes reversed                                     :  not reversed
        //    0    1     2     3     4     5     6      7     :     0     1     2     3    4     5      6      7
        //  len1  len3  atom1 atom3  a1    a2   a1-a2  a2-a3      len1  len3  atom1 atom3   a1    a3  a1-a2  a2-a3
    (rev ? [ h[2], h[0], h[5], h[3], h[8], h[6], h[10], h[9] ] : [ h[0], h[2], h[3], h[5], h[6], h[8],  h[9], h[10] ]);
    // h[1] = angle2 for both cases


/*
//
// generate the male or female interior cylinders of a rotating bond
//
*/
module joinUnit(cOuterLen, cOuterRad, cInnerLen, cInnerRad, male=false) {
    if (male) {
        colorSelect(bondColor) rotate([0,0,oRot]) {
            cylinder(h=cInnerLen, r=cInnerRad, center=false, $fn=8);
            cylinder(h=cOuterLen, r=cOuterRad, center=false, $fn=8);
        }
    } else {
        colorSelect(bondColor) rotate([0,0,fnRot]) {
            cylinder(h=cInnerLen, r=cInnerRad, center=false, $fn=30);
            cylinder(h=cOuterLen, r=cOuterRad, center=false, $fn=30);
        }
    }
}

/*
//
// create a rotatable bond
//
// supportSel : 0 for no support, 1 or 2 for support on top or bottom (needed
// for reversed hedra)
//
*/
module joiner(bondlen, scal, male=0, ver=0, supportSelect=0) {  // ver = differentiate joiner part lengths to guide assembly, but not used
    lenfac = bondLenFac;
    jClr = clearance+0.05;

    cOuterRad = (jBondRadius * scal) - (2*wall + (male ? jClr/2 : -jClr/2));
    cInnerRad = cOuterRad - joinerStep;  // m/f jClr already in cOuterRad;  - (male ? 0 : -0*jClr/2);

    hArmLen = (bondlen * lenfac);
    lenClr = 0.6*jClr;  // length clearance applied to male and female both, so effective clearance is 2x this value
    cOuterLen = hArmLen * lenfac + (ver ? 0.5 : - 0.5) - (wall+ (male ? lenClr*2 : -lenClr*2  ));

    joinerOffset = (hArmLen * (1 - lenfac)) + (male ? lenClr : -lenClr) - (ver ? 1 : 0);

    i=supportSelect-1;
    oside = cOuterRad*octSide;
    wid = oside+2*wall+4*jClr+1;

    if (male) {
        rotate([0,180,0])
        translate([0,0,-(bondlen-joinerOffset)]) {
            difference() {
                joinUnit(cOuterLen, cOuterRad, bondlen, cInnerRad, male=true);
                if (supportSelect) {
                    rotate([0,0,i*180]) {
                        translate([0,(cOuterRad*apmFac)-0.5*layerHeight,cOuterLen/2]) {
                            cube([oside+2*shim,layerHeight+shim,cOuterLen+2*shim],center=true);
                        }
                    }
                }
            }
            if (supportSelect) {
                rotate([0,0,i*180]) {
                    translate([0,(cOuterRad*apmFac)-0.5*layerHeight,cOuterLen/2]) {
                        for (j=[0:1]) {
                            rotate([0,(j?60:-60),0])
                                colorSelect(bondColor) cube([wid,layerHeight,2*nozzleDiameter],center=true);
                        }
                    }
                }
            }
        }
    } else {
        translate([0,0,joinerOffset]) {
            joinUnit(cOuterLen, cOuterRad, bondlen, cInnerRad);
            if (supportSelect) {  // extra gap top and bottom because filament sags
                supHeight = max(5*layerHeight,2*(cOuterRad-cOuterRad*apmFac));  // double because center=true below
                for(j=[0:1]) {
                    rotate([0,0,j*180]) {
                        translate([0,(cOuterRad*apmFac),cOuterLen/2]) {
                            colorSelect(bondColor) cube([oside+2*shim,supHeight+shim,cOuterLen+2*shim],center=true);
                        }
                    }
                }
            }
        }
    }
}


/*
//
// create bond with different options (regular, skinny, h-bond atom, rotatable
// male or female
//
//  parameters:
//  bl : bond length
//  br : bond radius
//  scal : protein_scale
//  key : option symbols defined below
//  atm : atomic element symbol, used for color and radius by atom() routine above
//  ver : make rotatable bonds slightly different based on value; currently unused
//  supporSel : enable print-in-place support for rotatable bonds
//
*/

// option symbols - these names generated in BioPython code so avoid changing without thought
StdBond = 1;
FemaleJoinBond = 2;
MaleJoinBond = 3;
SkinnyBond = 4;        // Calpha - Cbeta bond cylinder needs to be skinny for clearance with rotating bonds
HBond = 5;             // make room inside atom/bond to insert magnet to appropriate depth

module bond(bl, br, scal, key, atm, ver, supportSel=0) {

    br = (key == FemaleJoinBond ? jBondRadius * scal : br)  * (key == SkinnyBond ? 0.65 : 1);   // bond radius smaller for skinnyBond
    bl = (key == FemaleJoinBond ? bl * bondLenFac : bl);  // make female joiner shorter
    if (key == MaleJoinBond) { // male join is direct solid, others need difference()
        joiner(bl, scal, male = true, ver = ver, supportSelect=supportSel);
    } else {  // regular bond / skinny / h-bond / female join
        bhblen = bl +(hblen/2 * scal);
        rotate([0,0,fnRot]) {
            difference() {
                union() {
                    colorSelect(bondColor) cylinder(h=bl,r=br,center=false);
                    if (key == HBond) {  // make extension collar for h-bond magnet
                        rotate([0,0,oRot-fnRot]) 
                            colorSelect(bondColor) cylinder(h=bhblen-1,r=(magR + clearance +wall),center=false, $fn=8);
                    }
                }
                atom(atm,scal,-clearance);  // remove overlap with atom to clear area for female join
                if (key == HBond) {     // make space to insert magnet inside bond cylinder
                    translate([0,0,(bhblen-magL)-pClearance])
                        cylinder(h=magL+pClearance+shim, r=magR+pClearance, center=false, $fn=8);
                }
            }
        }
    }
}

/*
//
// Generate a 'hedron', one plane of 3 points, consisting of 3 atoms joined by
//   two bonds.
//   Defined as bond length - bond angle - bond length
//
// In some cases the sequence of atoms in the h[] array is reversed (rev flag),
// as detailed in the comments.
//
// other parameters:
//
// h = hedron array data according to rev flag:
//   yes reversed                                     :  not reversed
//    0    1     2     3     4     5     6      7     :     0     1     2     3    4     5      6      7
//  len1  len3  atom1 atom3  a1    a2   a1-a2  a2-a3      len1  len3  atom1 atom3   a1    a3  a1-a2  a2-a3
//
// split: chop half of the hedron - to selectively print parts of a rotating
//   bond to be glued together.  top or bottom half selected by global caTop
//   (C-alpha top) variable, undef by default so bottom half.
//
// supporSel: enable support structure inside rotatable bond to print in place.
//  Please note the bond needs to be exactly parallel to the buildplate and the
//  layerHeight global variable above needs to be set correctly for the
//  structure to be correctly created by your slicer software.
//
 */

module hedron(h,rev=0,scal,split=0, supportSel) {

    newh = hFlip(h, rev);  // make a consistent hedron array regardless of rev flag

    bondRad = bondRadius * scal;
    difference() {
        union(){
            if (h[7]) {
                // central atom at 0,0,0
                atom(h[4],scal);
            }

            if (newh[5] && newh[7] != FemaleJoinBond) {  // not female join
                // comments for non-reversed case
                // atom 3 is len3 up on +z
                translate([0,0,newh[1]])
                    difference() {
                        atom(newh[3],scal * (newh[7] == SkinnyBond ? 0.7 : 1));  // if skinny bond make atom (C-beta) same diameter as bond
                        if (newh[7] == HBond) {  // make room for hbond magnet through atom - this branch not used for backbone N,O
                            translate([0,0,scal*hblen/2-magL-pClearance])
                                cylinder(h=magL+pClearance,r=magR+pClearance,$fn=8);
                        }
                    }
            }

            if (newh[7]) {
                // atom 2 - atom 3 bond from origin up +z distance len3
                bond(newh[1], bondRad, scal, newh[7], h[4], ver=1, supportSel=supportSel);
            }
            rotate([0, h[1], 0]) {                        // rotate following elements by angle2 about Y
                if (newh[6]) {
                    bond(newh[0], bondRad, scal, newh[6], h[4], ver=1, supportSel=supportSel);  // h[4] is center atom (atom 2)
                }
                if (newh[4] && newh[6] != FemaleJoinBond) {   // if draw atom 2 and atom1-atom2 not joiner
                    translate([0,0,newh[0]]) {
                        difference() {
                            atom(newh[2],scal * (newh[6] == SkinnyBond ? 0.7 : 1));  // put atom1 sphere len1 away on Z
                            if (newh[6] == HBond) {  // make room for hbond magnet through atom
                                translate([0,0,scal*hblen/2-magL-pClearance])
                                    cylinder(h=magL+pClearance,r=magR+pClearance,$fn=8);
                            }
                        }
                    }
                }
            }
        }

        if (split) {
            // top / bottom half cutter
            thick = 2*bondRadius * scal;
            Zdim = newh[0];
            Xdim = newh[1];

            cside = 7* defaultAtomRadius * atomScale * scal / 12 + (caTop ? pClearance : -pClearance);
            difference() {
                translate([-Xdim,((rev || caTop) ? 0 : -thick),-Zdim]) {
                    colorSelect(bondColor) cube([2*Xdim,thick,2*Zdim]);
                }
                if (!caTop) {
                    rotate([0,(rev ? h[1] : 0),0])
                    rotate([45,0,0])
                    cube([cside, cside, cside],center=true);
                }
            }
            if (caTop) {
                //translate([tx+cside,0,tx+cside])
                    rotate([0,(rev ? h[1] : 0),0])
                        rotate([45,0,0])
                        colorSelect(bondColor) cube([cside, cside, cside], center=true);
            }
        }

        if (newh[7] == FemaleJoinBond) {  // female join
            joiner(newh[1], scal, male=false, ver=1, supportSelect=supportSel);
        }

        if (newh[6] == FemaleJoinBond) {  // female join
            rotate([0, h[1], 0]) {                        // rotate following elements by angle2 about Y
            joiner(newh[0], scal, male=false, ver=1, supportSelect=supportSel);
            translate([0,0,newh[0]])
                atom(newh[2],scal+0.5,clearance);  // clearance for atom against join outer cylinder
            }
        }

        if (newh[7] == FemaleJoinBond || newh[6] == FemaleJoinBond) {  // female join both hedron arms
            translate([0,0,newh[1]]) atom(newh[3],scal+0.5,clearance);  // clearance for atom against join outer cylinder
        }
    }
}

/*
//
// Hook to call custom routines for specific hedra.
//
// Residue is h[h_residue]
// Sequence position is h[h_seqpos]
//
*/
module hedronDispatch(h,rev=0,scal) {
    // default action is just to pass to hedron()

    hedron(h, rev, scal, 0, (support ? 1 : 0));

    /*
    // Some examples for special handling for specific hedra below:
    // note use of h_seqpos, h_residue, h_class for selecting hedra

    // bool flag caTop (for rotatable bond part) needs to be a global variable
    // so hedron() above can see it.

caBase1 = false;   // only make bottom of N_C-alpha_C hedron
caBase2 = false;   // same as caBase1 but for case of reversed hedron (for testing, should be identical to caBase1 result)
amideOnly = false; // make only the first amide

    if (caTop) {
        // these examples select a specific sequence position (h[h_seqpos] == n)
        if (h[h_seqpos] == 1) {
            if (h[h_class] == "NCAC") {
                hedron(h, rev, scal, 1);
            } else if (h[h_class] == "CBCAC") {
                colorSelect("yellow") {  // ca-cb
                    hedron(h, rev, scal);
                }
            }
        }
    } else if (caBase1) {
        if (h[h_seqpos] == 1 && (h[h_class] == "NCAC")) {
            hedron(h, rev, scal, true, (support ? 1 : 0));
        }
    } else if (caBase2) {
        if (h[h_seqpos] == 5 && (h[h_class] == "NCAC")) {
            hedron(h, rev, scal, true, (support ? 1 : 0));
        }
    } else if (amideOnly) {
        if (h[h_seqpos] == 1) {
            if (h[h_class] == "CACN") {
                colorSelect("darkgray") {
                    hedron(h, rev, scal);
                }
            }  else if (h[h_class] == "CACO") {
                colorSelect("red") {   // c=o
                    hedron(h, rev, scal);
                }
            }  else if (h[h_class] == "CNCA") {
                colorSelect("cyan") {  // h=n
                    hedron(h, rev, scal);
                }
            }
        } else if ((h[h_seqpos] == 2) && (h[h_class] == "HNCA")) {
            colorSelect("cyan") {  // h=n
                hedron(h, rev, scal);
            }
        }
       // actions above select out only a single hedron
    } else {
        // actions below will process hedra all but handle selected ones differently

        if (h[h_class] == "NCAC") {
            if (h[h_seqpos] == 1) {
                if (! CCap && NCap) {  // make split rotatable bond for terminal NH3
                    hedron(h, rev, scal, true, (support ? 1 : 0));
                }
            } else if (h[h_seqpos] == 5) {  // make split rotatable bond for terminal COOH
                hedron(h, rev, scal, true, (support ? 2 : 0));  // note supportSel = 2
            } else {
                hedron(h, rev, scal, 0, (support ? 2 : 0));
            }
        } else if (h[h_class] == "CBCAC") {
            colorSelect("yellow") {                     // ca-cb -- color yellow in OpenSCAD renderer
                if (h[h_seqpos] == 1 ) {         // don't make here for N-term
                } else if (h[h_seqpos] == 5 ) {  // don't make here for C-term
                } else {
                    hedron(h, rev, scal);       // otherwise do make here
                }
            }
        } else if (h[h_class] == "HNCA") {
            colorSelect("cyan") { // color h-n in OenSCAD renderer
                if (h[h_seqpos] == 1) {
                    if (NCap) {                      // only make at N term if variable NCap is true
                        hedron(h, rev, scal, 0, (support ? 1 : 0));
                    }
                } else {
                    hedron(h, rev, scal, 0, (support ? 1 : 0));
                }
            }
        } else if (h[h_residue] == "P") {
            colorSelect("darkgray")   // highlight Prolines in OpenSCAD renderer
                hedron(h, rev, scal);
        } else {
            echo("unrecognised hedron", h[h_class]);
            colorSelect("pink")
                hedron(h, rev, scal, 0, (support ? 1 : 0));
        }
    }
    */
}

/*
//
// Generate a hedron rotated to specific angle d
//
*/
module d2(d,hedra,scal)
{
    tz = (d[d_reversed] ? hedra[d[d_h2ndx]][2] : hedra[d[d_h2ndx]][0]);      // get h2 len1 depending on reversed
    rotate(d[d_dangle1]) {                                                   // 4. rotate h2 to specified dihedral angle
        translate([0,0,tz]) {                                               // 3. translate h2 h2:len1 up +z
            rotate([180, 0, 0]) {                                          // 2. rotate h2r about X so h2:a3 in +z and h2:a1 in -z
                hedronDispatch(hedra[d[d_h2ndx]],(!d[d_reversed]),scal);  // 1. reverse hedron 2 orientation = h2r
            }
        }
    }
}

/*
//
// Generate two hedra at specified dihedral angle d
//
*/
module dihedron(d,hedra,scal)
{
    if (d[d_h1new])
        hedronDispatch(hedra[d[d_h1ndx]],d[d_reversed],scal);                // reverse h1 if dihedral reversed
    if (d[d_h2new])
        d2(d,hedra,scal);
}

/*
//
// Generate a residue consisting of the set of dihedra in the parameter 'r',
//   referring to hedra the table specified in the parameter 'hedra'.
//
*/
module residue(r,hedra, scal)
{
    for (d = r) {
        multmatrix(d[d_dihedralTransform]) {
            dihedron(d, hedra, scal);
        }
    }
}

/*
//
// Generate a chain of residues, each positioned by a supplied
// rotation/translation matrix.
//
*/
module chain(protein)
{
    chnD = protein[p_chainData];
    c = chnD[c_residues];
    dihedra = chnD[c_dihedra];
    hedra = chnD[c_hedra];
    for (r = c) {
        multmatrix(r[r_resTransform]) {
            residue(dihedra[r[r_resNdx]],hedra, protein[p_proteinScale]);
        }
    }
}

/*
//
// OpenSCAD array indices to reference protein data - tied to BioPython code
//
*/

// protein base level
p_pdbid = 0;
p_proteinScale = 1;
p_chainData = 2;

// chain level data
c_chainID = 0;
c_dihedra = 1;
c_hedra = 2;
c_residues = 3;

// hedra definitions
h_len1 = 0;
h_angle2 = 1;
h_len3 = 2;
h_atom1class = 3;
h_atom2class = 4;
h_atom3class = 5;
h_atom1state = 6;
h_atom2state = 7;
h_atom3state = 8;
h_bond1state = 9;
h_bond2state = 10;
h_residue = 11;
h_seqpos = 12;  // residue sequence position for first atom in hedra
h_class = 13;

// dihedra specifications for each residue in sequence, dihedral array
d_dangle1 = 0;
d_h1ndx = 1;
d_h2ndx = 2;
d_reversed = 3;
d_h1new = 4;
d_h2new = 5;
d_dihedralTransform = 6;

// residueSet: world transform for each residue in sequence array
r_resNdx = 0;
r_resID = 1;
r_resTransform = 2;


// use single default atom radius for all atoms if tubes = true, else use
// covalent radii from literature
atomData = ( tubes ?
            [   ["Csb","green" , defaultAtomRadius], ["Cres","green" , defaultAtomRadius], ["Cdb","green" , defaultAtomRadius],
                ["Osb","red" , defaultAtomRadius], ["Ores","red" , defaultAtomRadius], ["Odb","red" , defaultAtomRadius],
                ["Nsb","blue" , defaultAtomRadius], ["Nres","blue" , defaultAtomRadius], ["Ndb","blue" , defaultAtomRadius],
                ["Hsb","gray" , defaultAtomRadius],
                ["Ssb","yellow" , defaultAtomRadius] ]
            :

// covalent radii from Heyrovska, Raji : 'Atomic Structures of all the Twenty
// Essential Amino Acids and a Tripeptide, with Bond Lengths as Sums of Atomic
// Covalent Radii'  https://arxiv.org/pdf/0804.2488.pdf

            [   ["Csb","green" , 0.77], ["Cres","green" , 0.72], ["Cdb","green" , 0.67],
                ["Osb","red" , 0.67], ["Ores","red" , 0.635], ["Odb","red" , 0.60],
                ["Nsb","blue" , 0.70], ["Nres","blue" , 0.66], ["Ndb","blue" , 0.62],
                ["Hsb","gray" , 0.37],
                ["Ssb","yellow" , 1.04] ]
    );


// optionally include protein array data here [ write_SCAD(includeCode=False) ], e.g.:
// include <1rtm.scad>;
// or paste below

protein = [ "0PDB", 10.0,  // ID, protein_scale
 [
   "A", // chain id
   [  // residue array of dihedra
     [ // 0 : (' ', 152, ' ') D backbone
      [ 152.34619, 0, 1, 0, 1, 1,     // 152_D_N:152_D_CA:152_D_C:152_D_O [ 152_D_N:152_D_CA:152_D_C -- 152_D_CA:152_D_C:152_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -26.52633, 0, 2, 0, 0, 1,     // 152_D_N:152_D_CA:152_D_C:153_I_N [ 152_D_N:152_D_CA:152_D_C -- 152_D_CA:152_D_C:153_I_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -178.91479, 2, 3, 0, 0, 1,     // 152_D_CA:152_D_C:153_I_N:153_I_CA [ 152_D_CA:152_D_C:153_I_N -- 152_D_C:153_I_N:153_I_CA ] 
        [ [ 0.3953943221940396, 0.4466089554021872, 0.8026230565659593, 0.0 ], [ -0.19736321807706675, 0.8947292556715507, -0.4006336471087982, 0.0 ], [ -0.8970569046203146, 0.0, 0.4419150482536432, 15.341225184500413 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.30307, 3, 4, 0, 0, 1,     // 152_D_C:153_I_N:153_I_CA:153_I_C [ 152_D_C:153_I_N:153_I_CA -- 153_I_N:153_I_CA:153_I_C ] 
        [ [ -0.893217834733885, -0.4390403323021219, 0.09698188658379553, 10.684238106259176 ], [ 0.4343559532038414, -0.8983067055080032, -0.06618133238111774, -5.333095335441121 ], [ 0.1161757531918433, -0.016989686618654205, 0.9930833524528079, 21.223844129540552 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 152, ' ') D sidechain

     ],
     [ // 1 : (' ', 153, ' ') I backbone
      [ -57.46720, 4, 5, 0, 0, 1,     // 153_I_N:153_I_CA:153_I_C:153_I_O [ 153_I_N:153_I_CA:153_I_C -- 153_I_CA:153_I_C:153_I_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 122.25606, 4, 6, 0, 0, 1,     // 153_I_N:153_I_CA:153_I_C:154_R_N [ 153_I_N:153_I_CA:153_I_C -- 153_I_CA:153_I_C:154_R_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.42864, 6, 7, 0, 0, 1,     // 153_I_CA:153_I_C:154_R_N:154_R_CA [ 153_I_CA:153_I_C:154_R_N -- 153_I_C:154_R_N:154_R_CA ] 
        [ [ -0.2334587914771561, -0.8456713751246809, -0.47993428506073305, 0.0 ], [ 0.36992308305330607, -0.5337039678461565, 0.76047155589779, 0.0 ], [ -0.8992518586616113, 0.0, 0.4374312456759775, 15.341356957211957 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -137.80914, 7, 8, 0, 0, 1,     // 153_I_C:154_R_N:154_R_CA:154_R_C [ 153_I_C:154_R_N:154_R_CA -- 154_R_N:154_R_CA:154_R_C ] 
        [ [ 0.5158038411620391, 0.8517552559701601, -0.09197489532308697, -6.396776893169154 ], [ -0.8458615188765365, 0.5233592507592328, 0.10302128677420617, 10.135901993465197 ], [ 0.13588483479168895, 0.024659249216905266, 0.9904177063752919, 21.171634501014026 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 153, ' ') I sidechain

     ],
     [ // 2 : (' ', 154, ' ') R backbone
      [ -43.09779, 8, 9, 0, 0, 1,     // 154_R_N:154_R_CA:154_R_C:154_R_O [ 154_R_N:154_R_CA:154_R_C -- 154_R_CA:154_R_C:154_R_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 136.40588, 8, 10, 0, 0, 1,     // 154_R_N:154_R_CA:154_R_C:155_Q_N [ 154_R_N:154_R_CA:154_R_C -- 154_R_CA:154_R_C:155_Q_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -178.78102, 10, 11, 0, 0, 1,     // 154_R_CA:154_R_C:155_Q_N:155_Q_CA [ 154_R_CA:154_R_C:155_Q_N -- 154_R_C:155_Q_N:155_Q_CA ] 
        [ [ -0.32849558138276275, -0.6895452437754972, -0.645459533819581, 0.0 ], [ 0.3127578561497455, -0.7242426090659055, 0.614536546211879, 0.0 ], [ -0.8912200493865787, 0.0, 0.45357118908875177, 15.235613145305647 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -68.26433, 11, 12, 0, 0, 1,     // 154_R_C:155_Q_N:155_Q_CA:155_Q_C [ 154_R_C:155_Q_N:155_Q_CA -- 155_Q_N:155_Q_CA:155_Q_C ] 
        [ [ 0.729672769851978, 0.682400893634347, -0.04366542457811507, -8.58586118025412 ], [ -0.6794037460923721, 0.7307322065318035, 0.0666407692993818, 8.174525589148482 ], [ 0.07738345257349256, -0.01895950168666127, 0.9968211166320663, 21.26898791294422 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 154, ' ') R sidechain

     ],
     [ // 3 : (' ', 155, ' ') Q backbone
      [ -53.46226, 12, 13, 0, 0, 1,     // 155_Q_N:155_Q_CA:155_Q_C:155_Q_O [ 155_Q_N:155_Q_CA:155_Q_C -- 155_Q_CA:155_Q_C:155_Q_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 126.97448, 12, 14, 0, 0, 1,     // 155_Q_N:155_Q_CA:155_Q_C:156_G_N [ 155_Q_N:155_Q_CA:155_Q_C -- 155_Q_CA:155_Q_C:156_G_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.29436, 14, 15, 0, 0, 1,     // 155_Q_CA:155_Q_C:156_G_N:156_G_CA [ 155_Q_CA:155_Q_C:156_G_N -- 155_Q_C:156_G_N:156_G_CA ] 
        [ [ -0.2674087959250598, -0.7989034420817724, -0.5387446761609862, 0.0 ], [ 0.35519245825716245, -0.6014593005598102, 0.7156011650125633, 0.0 ], [ -0.8957292299903714, 0.0, 0.4445999848637608, 15.203618227152155 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -66.89920, 15, 16, 0, 0, 1,     // 155_Q_C:156_G_N:156_G_CA:156_G_C [ 155_Q_C:156_G_N:156_G_CA -- 156_G_N:156_G_CA:156_G_C ] 
        [ [ 0.5880196130923234, 0.8065088062116518, -0.06145307251725233, -7.164170199181778 ], [ -0.8060114605584737, 0.5906206335242771, 0.03889463592550276, 9.515989238937417 ], [ 0.06766431901046698, 0.02666107196715544, 0.9973518572572129, 21.115862195509877 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 155, ' ') Q sidechain

     ],
     [ // 4 : (' ', 156, ' ') G backbone
      [ -28.05505, 16, 17, 0, 0, 1,     // 156_G_N:156_G_CA:156_G_C:156_G_O [ 156_G_N:156_G_CA:156_G_C -- 156_G_CA:156_G_C:156_G_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 150.97703, 16, 18, 0, 0, 1,     // 156_G_N:156_G_CA:156_G_C:157_P_N [ 156_G_N:156_G_CA:156_G_C -- 156_G_CA:156_G_C:157_P_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.96915, 18, 19, 0, 0, 1,     // 156_G_CA:156_G_C:157_P_N:157_P_CA [ 156_G_CA:156_G_C:157_P_N -- 156_G_C:157_P_N:157_P_CA ] 
        [ [ -0.420539632964679, -0.4851602206863669, -0.7666589707095907, 0.0 ], [ 0.23332937318093946, -0.8744252742593593, 0.4253678917683305, 0.0 ], [ -0.8767575609693498, 0.0, 0.480932614077146, 15.229125383188661 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -56.91524, 19, 20, 0, 0, 1,     // 156_G_C:157_P_N:157_P_CA:157_P_C [ 156_G_C:157_P_N:157_P_CA -- 157_P_N:157_P_CA:157_P_C ] 
        [ [ 0.8733003906991935, 0.4849337274313226, -0.04675155189081257, -10.30435079500439 ], [ -0.48421140036990096, 0.8745507744912406, 0.026462475103256427, 5.717196486536209 ], [ 0.053719152603658275, -0.0004720554360250758, 0.9985559722955996, 21.693144658437696 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 156, ' ') G sidechain

     ],
     [ // 5 : (' ', 157, ' ') P backbone
      [ 156.18449, 20, 21, 0, 0, 1,     // 157_P_N:157_P_CA:157_P_C:157_P_O [ 157_P_N:157_P_CA:157_P_C -- 157_P_CA:157_P_C:157_P_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -23.47457, 20, 22, 0, 0, 1,     // 157_P_N:157_P_CA:157_P_C:158_K_N [ 157_P_N:157_P_CA:157_P_C -- 157_P_CA:157_P_C:158_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.08250, 22, 23, 0, 0, 1,     // 157_P_CA:157_P_C:158_K_N:158_K_CA [ 157_P_CA:157_P_C:158_K_N -- 157_P_C:158_K_N:158_K_CA ] 
        [ [ 0.4177913877837107, 0.39834204080251473, 0.8165622908407225, 0.0 ], [ -0.1814404385113306, 0.91723694786532, -0.3546205701077257, 0.0 ], [ -0.8902413849999208, 0.0, 0.45548905193585376, 15.322276579842605 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -120.90687, 23, 24, 0, 0, 1,     // 157_P_C:158_K_N:158_K_CA:158_K_C [ 157_P_C:158_K_N:158_K_CA -- 158_K_N:158_K_CA:158_K_C ] 
        [ [ -0.9105043150474754, -0.40498095493371733, 0.0835004097049949, 10.866608735188324 ], [ 0.40466343056620835, -0.9142139876467463, -0.0214544343533519, -4.719202721010463 ], [ 0.08502587983852367, 0.014255207189277452, 0.996276763166579, 21.383812075747915 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 157, ' ') P sidechain

     ],
     [ // 6 : (' ', 158, ' ') K backbone
      [ -170.24700, 24, 25, 0, 0, 1,     // 158_K_N:158_K_CA:158_K_C:158_K_O [ 158_K_N:158_K_CA:158_K_C -- 158_K_CA:158_K_C:158_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  10.25075, 24, 26, 0, 0, 1,     // 158_K_N:158_K_CA:158_K_C:159_E_N [ 158_K_N:158_K_CA:158_K_C -- 158_K_CA:158_K_C:159_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.89412, 26, 27, 0, 0, 1,     // 158_K_CA:158_K_C:159_E_N:159_E_CA [ 158_K_CA:158_K_C:159_E_N -- 158_K_C:159_E_N:159_E_CA ] 
        [ [ 0.43652082000886955, -0.17795646095121045, 0.8819189711671388, 0.0 ], [ 0.0789417397281942, 0.9840383620599963, 0.15948888275665404, 0.0 ], [ -0.896224176993385, 0.0, 0.44360142534997543, 15.299185513032217 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -83.89681, 27, 28, 0, 0, 1,     // 158_K_C:159_E_N:159_E_CA:159_E_C [ 158_K_C:159_E_N:159_E_CA -- 159_E_N:159_E_CA:159_E_C ] 
        [ [ -0.979299092697733, 0.17714952142926718, 0.09793535673483163, 11.738264379906177 ], [ -0.17610256890290413, -0.9841825564181145, 0.019302353952374028, 2.1227830817335542 ], [ 0.09980567252012798, 0.0016561098050642048, 0.9950055703528089, 21.203481690231005 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 158, ' ') K sidechain

     ],
     [ // 7 : (' ', 159, ' ') E backbone
      [ -51.50291, 28, 29, 0, 0, 1,     // 159_E_N:159_E_CA:159_E_C:159_E_O [ 159_E_N:159_E_CA:159_E_C -- 159_E_CA:159_E_C:159_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 127.73064, 28, 30, 0, 0, 1,     // 159_E_N:159_E_CA:159_E_C:160_P_N [ 159_E_N:159_E_CA:159_E_C -- 159_E_CA:159_E_C:160_P_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.69246, 30, 31, 0, 0, 1,     // 159_E_CA:159_E_C:160_P_N:160_P_CA [ 159_E_CA:159_E_C:160_P_N -- 159_E_C:160_P_N:160_P_CA ] 
        [ [ -0.28392308974959723, -0.7908963887715672, -0.5420983133483593, 0.0 ], [ 0.36694781676852717, -0.6119500814838527, 0.7006185820692435, 0.0 ], [ -0.8858538134906083, 0.0, 0.46396446105725214, 15.263489732410433 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -58.10197, 31, 32, 0, 0, 1,     // 159_E_C:160_P_N:160_P_CA:160_P_C [ 159_E_C:160_P_N:160_P_CA -- 160_P_N:160_P_CA:160_P_C ] 
        [ [ 0.6085676655181514, 0.7924089658415642, -0.04163444895384121, -7.302255986873195 ], [ -0.7910690562332009, 0.6099716568417485, 0.046306869037380585, 9.437580065190266 ], [ 0.06208981201533922, 0.004754861053171072, 0.9980592399954348, 21.513255097894145 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 159, ' ') E sidechain

     ],
     [ // 8 : (' ', 160, ' ') P backbone
      [ -37.83347, 32, 33, 0, 0, 1,     // 160_P_N:160_P_CA:160_P_C:160_P_O [ 160_P_N:160_P_CA:160_P_C -- 160_P_CA:160_P_C:160_P_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 140.92863, 32, 34, 0, 0, 1,     // 160_P_N:160_P_CA:160_P_C:161_F_N [ 160_P_N:160_P_CA:160_P_C -- 160_P_CA:160_P_C:161_F_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.09705, 34, 35, 0, 0, 1,     // 160_P_CA:160_P_C:161_F_N:161_F_CA [ 160_P_CA:160_P_C:161_F_N -- 160_P_C:161_F_N:161_F_CA ] 
        [ [ -0.35042173307546914, -0.6302879579642553, -0.6927782466515774, 0.0 ], [ 0.28448939676487756, -0.7763614429151213, 0.5624310562932198, 0.0 ], [ -0.892339841157359, 0.0, 0.4513641632687062, 15.288478941330302 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -54.40090, 35, 36, 0, 0, 1,     // 160_P_C:161_F_N:161_F_CA:161_F_C [ 160_P_C:161_F_N:161_F_CA -- 161_F_N:161_F_CA:161_F_C ] 
        [ [ 0.7787794451499263, 0.6246874616530158, -0.057167744974635754, -9.214516419149113 ], [ -0.6215917362330381, 0.7807482497299468, 0.06368582252292702, 7.480792342863112 ], [ 0.08441735164508657, -0.014062211673010397, 0.9963312526183715, 21.291990907182708 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 160, ' ') P sidechain

     ],
     [ // 9 : (' ', 161, ' ') F backbone
      [ 137.08994, 36, 37, 0, 0, 1,     // 161_F_N:161_F_CA:161_F_C:161_F_O [ 161_F_N:161_F_CA:161_F_C -- 161_F_CA:161_F_C:161_F_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.34191, 36, 38, 0, 0, 1,     // 161_F_N:161_F_CA:161_F_C:162_R_N [ 161_F_N:161_F_CA:161_F_C -- 161_F_CA:161_F_C:162_R_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.94706, 38, 39, 0, 0, 1,     // 161_F_CA:161_F_C:162_R_N:162_R_CA [ 161_F_CA:161_F_C:162_R_N -- 161_F_C:162_R_N:162_R_CA ] 
        [ [ 0.33463532468893203, 0.6735533442404578, 0.659048626379593, 0.0 ], [ -0.3049424524895067, 0.7391386151883119, -0.6005699028482775, 0.0 ], [ -0.8916441555575957, 0.0, 0.45273689915886245, 15.318934157581468 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.27735, 39, 40, 0, 0, 1,     // 161_F_C:162_R_N:162_R_CA:162_R_C [ 161_F_C:162_R_N:162_R_CA -- 162_R_N:162_R_CA:162_R_C ] 
        [ [ -0.7360439016244039, -0.6738622425354777, 0.064413142809382, 8.77456846752067 ], [ 0.6713933844706982, -0.7388565485828991, -0.05763613367625637, -7.995983181124526 ], [ 0.08643088666967405, 0.0008238332495679596, 0.9962575084426069, 21.346669836609124 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 161, ' ') F sidechain

     ],
     [ // 10 : (' ', 162, ' ') R backbone
      [ 141.15041, 40, 41, 0, 0, 1,     // 162_R_N:162_R_CA:162_R_C:162_R_O [ 162_R_N:162_R_CA:162_R_C -- 162_R_CA:162_R_C:162_R_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -37.94600, 40, 42, 0, 0, 1,     // 162_R_N:162_R_CA:162_R_C:163_D_N [ 162_R_N:162_R_CA:162_R_C -- 162_R_CA:162_R_C:163_D_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.14042, 42, 43, 0, 0, 1,     // 162_R_CA:162_R_C:163_D_N:163_D_CA [ 162_R_CA:162_R_C:163_D_N -- 162_R_C:163_D_N:163_D_CA ] 
        [ [ 0.3581713882601918, 0.6149185830928227, 0.7025584621929213, 0.0 ], [ -0.27929098097530175, 0.7885906011151258, -0.5478333795769825, 0.0 ], [ -0.8909039255596647, 0.0, 0.45419180466228076, 15.35243598600596 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -66.49005, 43, 44, 0, 0, 1,     // 162_R_C:163_D_N:163_D_CA:163_D_C [ 162_R_C:163_D_N:163_D_CA -- 163_D_N:163_D_CA:163_D_C ] 
        [ [ -0.7807066593954455, -0.6202226246038098, 0.07629553004708436, 9.355965670589372 ], [ 0.6188495603063752, -0.7843119702433422, -0.04335844833031603, -7.295492928128378 ], [ 0.08673138811416614, 0.013365225870448937, 0.9961420767406735, 21.4009047768017 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 162, ' ') R sidechain

     ],
     [ // 11 : (' ', 163, ' ') D backbone
      [ 139.48940, 44, 45, 0, 0, 1,     // 163_D_N:163_D_CA:163_D_C:163_D_O [ 163_D_N:163_D_CA:163_D_C -- 163_D_CA:163_D_C:163_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -40.76307, 44, 46, 0, 0, 1,     // 163_D_N:163_D_CA:163_D_C:164_Y_N [ 163_D_N:163_D_CA:163_D_C -- 163_D_CA:163_D_C:164_Y_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.15443, 46, 47, 0, 0, 1,     // 163_D_CA:163_D_C:164_Y_N:164_Y_CA [ 163_D_CA:163_D_C:164_Y_N -- 163_D_C:164_Y_N:164_Y_CA ] 
        [ [ 0.3265453718701161, 0.65293253283591, 0.6834085364368014, 0.0 ], [ -0.2814993031641122, 0.7574160729509793, -0.5891341398642893, 0.0 ], [ -0.9022894560108341, 0.0, 0.431130766208668, 15.296415621096306 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -65.51139, 47, 48, 0, 0, 1,     // 163_D_C:164_Y_N:164_Y_CA:164_Y_C [ 163_D_C:164_Y_N:164_Y_CA -- 164_Y_N:164_Y_CA:164_Y_C ] 
        [ [ -0.7397587008870111, -0.6631104493421288, 0.114199809260969, 9.073274161545353 ], [ 0.6607763651167197, -0.7479572950686361, -0.06272543389044812, -7.821640035087882 ], [ 0.1270104710844636, 0.0290588493731138, 0.9914756293061435, 21.02032369763831 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 163, ' ') D sidechain

     ],
     [ // 12 : (' ', 164, ' ') Y backbone
      [ 135.21580, 48, 49, 0, 0, 1,     // 164_Y_N:164_Y_CA:164_Y_C:164_Y_O [ 164_Y_N:164_Y_CA:164_Y_C -- 164_Y_CA:164_Y_C:164_Y_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.39216, 48, 50, 0, 0, 1,     // 164_Y_N:164_Y_CA:164_Y_C:165_V_N [ 164_Y_N:164_Y_CA:164_Y_C -- 164_Y_CA:164_Y_C:165_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.15049, 50, 51, 0, 0, 1,     // 164_Y_CA:164_Y_C:165_V_N:165_V_CA [ 164_Y_CA:164_Y_C:165_V_N -- 164_Y_C:165_V_N:165_V_CA ] 
        [ [ 0.328391956454584, 0.6742013008154498, 0.6615219791622088, 0.0 ], [ -0.29978064294456025, 0.7385476328435122, -0.6038865457492746, 0.0 ], [ -0.8957065864733141, 0.0, 0.4446456015168972, 15.315436383197017 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -57.96507, 51, 52, 0, 0, 1,     // 164_Y_C:165_V_N:165_V_CA:165_V_C [ 164_Y_C:165_V_N:165_V_CA -- 165_V_N:165_V_CA:165_V_C ] 
        [ [ -0.7292309632823879, -0.6789960266692687, 0.08477380466621873, 8.81442163827133 ], [ 0.6764323962918264, -0.7340218235652726, -0.060424959882563516, -8.046460742928888 ], [ 0.10325413036357467, 0.013279996131684103, 0.9945663508613211, 21.240097898353838 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 164, ' ') Y sidechain

     ],
     [ // 13 : (' ', 165, ' ') V backbone
      [ 134.41814, 52, 53, 0, 0, 1,     // 165_V_N:165_V_CA:165_V_C:165_V_O [ 165_V_N:165_V_CA:165_V_C -- 165_V_CA:165_V_C:165_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -43.96463, 52, 54, 0, 0, 1,     // 165_V_N:165_V_CA:165_V_C:166_D_N [ 165_V_N:165_V_CA:165_V_C -- 165_V_CA:165_V_C:166_D_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.55680, 54, 55, 0, 0, 1,     // 165_V_CA:165_V_C:166_D_N:166_D_CA [ 165_V_CA:165_V_C:166_D_N -- 165_V_C:166_D_N:166_D_CA ] 
        [ [ 0.32222651313473505, 0.6942142189837959, 0.6436122220667881, 0.0 ], [ -0.3107863756705995, 0.7197684476001419, -0.620761798548215, 0.0 ], [ -0.8941934370876155, 0.0, 0.44768079819156487, 15.25645574170156 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.42939, 55, 56, 0, 0, 1,     // 165_V_C:166_D_N:166_D_CA:166_D_C [ 165_V_C:166_D_N:166_D_CA -- 166_D_N:166_D_CA:166_D_C ] 
        [ [ -0.7147685126910096, -0.6966859590953156, 0.06111176371961651, 8.567284888929498 ], [ 0.6949425410079196, -0.7173428970047702, -0.04973965032305363, -8.26311712858514 ], [ 0.07849100561808296, 0.00691682848160596, 0.9968908262797979, 21.215647947676302 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 165, ' ') V sidechain

     ],
     [ // 14 : (' ', 166, ' ') D backbone
      [ 136.78626, 56, 57, 0, 0, 1,     // 166_D_N:166_D_CA:166_D_C:166_D_O [ 166_D_N:166_D_CA:166_D_C -- 166_D_CA:166_D_C:166_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.20923, 56, 58, 0, 0, 1,     // 166_D_N:166_D_CA:166_D_C:167_R_N [ 166_D_N:166_D_CA:166_D_C -- 166_D_CA:166_D_C:167_R_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.81223, 58, 59, 0, 0, 1,     // 166_D_CA:166_D_C:167_R_N:167_R_CA [ 166_D_CA:166_D_C:167_R_N -- 166_D_C:167_R_N:167_R_CA ] 
        [ [ 0.342531180411564, 0.6718399193192495, 0.6567370198602823, 0.0 ], [ -0.31068887033490916, 0.7406963769380166, -0.5956855739741509, 0.0 ], [ -0.8866480791700158, 0.0, 0.4624447898983417, 15.30165526918609 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.10754, 59, 60, 0, 0, 1,     // 166_D_C:167_R_N:167_R_CA:167_R_C [ 166_D_C:167_R_N:167_R_CA -- 167_R_N:167_R_CA:167_R_C ] 
        [ [ -0.7325157210061076, -0.6787959062668429, 0.05154450615029373, 8.733652505041649 ], [ 0.6786663399834493, -0.7340969393712129, -0.022664566602677825, -7.921756575353897 ], [ 0.0532232792335412, 0.01837936998902274, 0.9984134821337476, 21.451502126606055 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 166, ' ') D sidechain

     ],
     [ // 15 : (' ', 167, ' ') R backbone
      [ 133.90448, 60, 61, 0, 0, 1,     // 167_R_N:167_R_CA:167_R_C:167_R_O [ 167_R_N:167_R_CA:167_R_C -- 167_R_CA:167_R_C:167_R_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -44.23739, 60, 62, 0, 0, 1,     // 167_R_N:167_R_CA:167_R_C:168_F_N [ 167_R_N:167_R_CA:167_R_C -- 167_R_CA:167_R_C:168_F_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.61711, 62, 63, 0, 0, 1,     // 167_R_CA:167_R_C:168_F_N:168_F_CA [ 167_R_CA:167_R_C:168_F_N -- 167_R_C:168_F_N:168_F_CA ] 
        [ [ 0.3083991107566697, 0.6976327500347652, 0.6466827155285861, 0.0 ], [ -0.3002968176342759, 0.7164555436863692, -0.6296930007867291, 0.0 ], [ -0.9026138763631001, 0.0, 0.4304511472824508, 15.259098274125872 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -57.26418, 63, 64, 0, 0, 1,     // 167_R_C:168_F_N:168_F_CA:168_F_C [ 167_R_C:168_F_N:168_F_CA -- 168_F_N:168_F_CA:168_F_C ] 
        [ [ -0.7068134733524439, -0.6996781059452298, 0.10423656723218394, 8.590228735823926 ], [ 0.6933985789797265, -0.714432757123611, -0.09372430964083717, -8.364545363925453 ], [ 0.14004686557131502, 0.0060318827617209915, 0.9901265029449514, 20.977008658867913 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 167, ' ') R sidechain

     ],
     [ // 16 : (' ', 168, ' ') F backbone
      [ 130.25011, 64, 65, 0, 0, 1,     // 168_F_N:168_F_CA:168_F_C:168_F_O [ 168_F_N:168_F_CA:168_F_C -- 168_F_CA:168_F_C:168_F_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -49.38859, 64, 66, 0, 0, 1,     // 168_F_N:168_F_CA:168_F_C:169_Y_N [ 168_F_N:168_F_CA:168_F_C -- 168_F_CA:168_F_C:169_Y_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.83913, 66, 67, 0, 0, 1,     // 168_F_CA:168_F_C:169_Y_N:169_Y_CA [ 168_F_CA:168_F_C:169_Y_N -- 168_F_C:169_Y_N:169_Y_CA ] 
        [ [ 0.2887278054229685, 0.7591416709356426, 0.5833868166359923, 0.0 ], [ -0.3367287502929609, 0.6509254361668781, -0.6803747682543057, 0.0 ], [ -0.8962421565078141, 0.0, 0.4435650988279204, 15.253761061450023 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -70.00653, 67, 68, 0, 0, 1,     // 168_F_C:169_Y_N:169_Y_CA:169_Y_C [ 168_F_C:169_Y_N:169_Y_CA -- 169_Y_N:169_Y_CA:169_Y_C ] 
        [ [ -0.6491080104643883, -0.758328012899291, 0.05997845949319022, 7.7636705742434575 ], [ 0.7547504848871105, -0.6518683092579514, -0.07361666216184612, -9.05437938795509 ], [ 0.09492363411519364, -0.002516393727637716, 0.9954813767464323, 21.156694085184135 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 168, ' ') F sidechain

     ],
     [ // 17 : (' ', 169, ' ') Y backbone
      [ 156.32812, 68, 69, 0, 0, 1,     // 169_Y_N:169_Y_CA:169_Y_C:169_Y_O [ 169_Y_N:169_Y_CA:169_Y_C -- 169_Y_CA:169_Y_C:169_Y_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -23.65104, 68, 70, 0, 0, 1,     // 169_Y_N:169_Y_CA:169_Y_C:170_K_N [ 169_Y_N:169_Y_CA:169_Y_C -- 169_Y_CA:169_Y_C:170_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.36395, 70, 71, 0, 0, 1,     // 169_Y_CA:169_Y_C:170_K_N:170_K_CA [ 169_Y_CA:169_Y_C:170_K_N -- 169_Y_C:170_K_N:170_K_CA ] 
        [ [ 0.4053639414564411, 0.401165192441122, 0.821429585138359, 0.0 ], [ -0.17752935285356533, 0.916005725076911, -0.3597455218647158, 0.0 ], [ -0.8967515842429793, 0.0, 0.4425342881152903, 15.192994371301106 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -67.59483, 71, 72, 0, 0, 1,     // 169_Y_C:170_K_N:170_K_CA:170_K_C [ 169_Y_C:170_K_N:170_K_CA -- 170_K_N:170_K_CA:170_K_C ] 
        [ [ -0.9061242248461814, -0.41257500681870396, 0.09338497146464976, 10.894327231625649 ], [ 0.4131204992301756, -0.9105637619275231, -0.014320914084216587, -4.771176381047949 ], [ 0.0909414221502916, 0.02560271885842229, 0.9955270757364357, 21.062168592887517 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 169, ' ') Y sidechain

     ],
     [ // 18 : (' ', 170, ' ') K backbone
      [ 133.75064, 72, 73, 0, 0, 1,     // 170_K_N:170_K_CA:170_K_C:170_K_O [ 170_K_N:170_K_CA:170_K_C -- 170_K_CA:170_K_C:170_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.85313, 72, 74, 0, 0, 1,     // 170_K_N:170_K_CA:170_K_C:171_T_N [ 170_K_N:170_K_CA:170_K_C -- 170_K_CA:170_K_C:171_T_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.61740, 74, 75, 0, 0, 1,     // 170_K_CA:170_K_C:171_T_N:171_T_CA [ 170_K_CA:170_K_C:171_T_N -- 170_K_C:171_T_N:171_T_CA ] 
        [ [ 0.3242227328495406, 0.6801213600018896, 0.6575063156903935, 0.0 ], [ -0.3007924489822243, 0.7330995400825048, -0.6099909564625495, 0.0 ], [ -0.8968854565321324, 0.0, 0.44226290581638045, 15.301235140151578 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -69.66573, 75, 76, 0, 0, 1,     // 170_K_C:171_T_N:171_T_CA:171_T_C [ 170_K_C:171_T_N:171_T_CA -- 171_T_N:171_T_CA:171_T_C ] 
        [ [ -0.7205779046120724, -0.6877463836573261, 0.08816118846264082, 8.75905934076906 ], [ 0.6859728282895816, -0.7256284157045977, -0.053895094122079434, -8.126077054298339 ], [ 0.10103841959011185, 0.02164056580373014, 0.9946471352587444, 21.192900275259923 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 170, ' ') K sidechain

     ],
     [ // 19 : (' ', 171, ' ') T backbone
      [ 144.19621, 76, 77, 0, 0, 1,     // 171_T_N:171_T_CA:171_T_C:171_T_O [ 171_T_N:171_T_CA:171_T_C -- 171_T_CA:171_T_C:171_T_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -34.68578, 76, 78, 0, 0, 1,     // 171_T_N:171_T_CA:171_T_C:172_L_N [ 171_T_N:171_T_CA:171_T_C -- 171_T_CA:171_T_C:172_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.85082, 78, 79, 0, 0, 1,     // 171_T_CA:171_T_C:172_L_N:172_L_CA [ 171_T_CA:171_T_C:172_L_N -- 171_T_C:172_L_N:172_L_CA ] 
        [ [ 0.36251521346776167, 0.5690755080175057, 0.7380621831384134, 0.0 ], [ -0.25088437876134706, 0.8222852705567684, -0.5107875902162416, 0.0 ], [ -0.8975743693410346, 0.0, 0.44086307568455296, 15.291600823145714 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -68.82921, 79, 80, 0, 0, 1,     // 171_T_C:172_L_N:172_L_CA:172_L_C [ 171_T_C:172_L_N:172_L_CA -- 172_L_N:172_L_CA:172_L_C ] 
        [ [ -0.8177926159762119, -0.5700174435125966, 0.07934325016123177, 9.80398737219825 ], [ 0.5676257883836554, -0.8216292671550862, -0.0522140949938748, -6.785004297417676 ], [ 0.09495368142740793, 0.002336973589819486, 0.9954789485156493, 21.147768671362464 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 171, ' ') T sidechain

     ],
     [ // 20 : (' ', 172, ' ') L backbone
      [ 140.23015, 80, 81, 0, 0, 1,     // 172_L_N:172_L_CA:172_L_C:172_L_O [ 172_L_N:172_L_CA:172_L_C -- 172_L_CA:172_L_C:172_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -37.90207, 80, 82, 0, 0, 1,     // 172_L_N:172_L_CA:172_L_C:173_R_N [ 172_L_N:172_L_CA:172_L_C -- 172_L_CA:172_L_C:173_R_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.35677, 82, 83, 0, 0, 1,     // 172_L_CA:172_L_C:173_R_N:173_R_CA [ 172_L_CA:172_L_C:173_R_N -- 172_L_C:173_R_N:173_R_CA ] 
        [ [ 0.3571966859809482, 0.6143137647565946, 0.7035830625838055, 0.0 ], [ -0.2780908008045787, 0.7890618470250474, -0.5477653768515209, 0.0 ], [ -0.8916703617548898, 0.0, 0.4526852835777898, 15.323498548437428 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.92717, 83, 84, 0, 0, 1,     // 172_L_C:173_R_N:173_R_CA:173_R_C [ 172_L_C:173_R_N:173_R_CA -- 173_R_N:173_R_CA:173_R_C ] 
        [ [ -0.7826290928345138, -0.6182850419166561, 0.07221710317608226, 9.379341323696774 ], [ 0.616789817400441, -0.7858901984752142, -0.044123883461444545, -7.3021633237253285 ], [ 0.0840358506838412, 0.010010138995398437, 0.9964124512053911, 21.358165971070015 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 172, ' ') L sidechain

     ],
     [ // 21 : (' ', 173, ' ') R backbone
      [ 146.81178, 84, 85, 0, 0, 1,     // 173_R_N:173_R_CA:173_R_C:173_R_O [ 173_R_N:173_R_CA:173_R_C -- 173_R_CA:173_R_C:173_R_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -33.96952, 84, 86, 0, 0, 1,     // 173_R_N:173_R_CA:173_R_C:174_A_N [ 173_R_N:173_R_CA:173_R_C -- 173_R_CA:173_R_C:174_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.62777, 86, 87, 0, 0, 1,     // 173_R_CA:173_R_C:174_A_N:174_A_CA [ 173_R_CA:173_R_C:174_A_N -- 173_R_C:174_A_N:174_A_CA ] 
        [ [ 0.36720993822580733, 0.5587518228822088, 0.7436082716686303, 0.0 ], [ -0.24740212762570793, 0.8293349145103013, -0.500994797198997, 0.0 ], [ -0.8966320586028985, 0.0, 0.44277641252163413, 15.210729581753663 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.66040, 87, 88, 0, 0, 1,     // 173_R_C:174_A_N:174_A_CA:174_A_C [ 173_R_C:174_A_N:174_A_CA -- 174_A_N:174_A_CA:174_A_C ] 
        [ [ -0.8264752599080624, -0.5563543891062817, 0.08607228637649263, 9.89890054101432 ], [ 0.5526218772803146, -0.8309247027618537, -0.06460030256357897, -6.66923413575003 ], [ 0.10746025084226124, -0.005825123372187198, 0.9941923166201885, 21.104961580576767 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 173, ' ') R sidechain

     ],
     [ // 22 : (' ', 174, ' ') A backbone
      [ 161.31868, 88, 89, 0, 0, 1,     // 174_A_N:174_A_CA:174_A_C:174_A_O [ 174_A_N:174_A_CA:174_A_C -- 174_A_CA:174_A_C:174_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -18.74441, 88, 90, 0, 0, 1,     // 174_A_N:174_A_CA:174_A_C:175_E_N [ 174_A_N:174_A_CA:174_A_C -- 174_A_CA:174_A_C:175_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.96989, 90, 91, 0, 0, 1,     // 174_A_CA:174_A_C:175_E_N:175_E_CA [ 174_A_CA:174_A_C:175_E_N -- 174_A_C:175_E_N:175_E_CA ] 
        [ [ 0.4225288545687902, 0.3213470244584229, 0.8474700330563341, 0.0 ], [ -0.1433832205470571, 0.9469615039016728, -0.2875850521043057, 0.0 ], [ -0.8949360977870655, 0.0, 0.4461943308443753, 15.24342077416869 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -104.16653, 91, 92, 0, 0, 1,     // 174_A_C:175_E_N:175_E_CA:175_E_C [ 174_A_C:175_E_N:175_E_CA -- 175_E_N:175_E_CA:175_E_C ] 
        [ [ -0.9425619272243125, -0.32156902627428974, 0.09039012494844333, 11.265641161147894 ], [ 0.32014786838424003, -0.9468860227898829, -0.030202685545903586, -3.8229434362804966 ], [ 0.09530139409377685, 0.0004703043297324789, 0.9954483528027057, 21.17479890361471 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 174, ' ') A sidechain

     ],
     [ // 23 : (' ', 175, ' ') E backbone
      [ -175.95026, 92, 93, 0, 0, 1,     // 175_E_N:175_E_CA:175_E_C:175_E_O [ 175_E_N:175_E_CA:175_E_C -- 175_E_CA:175_E_C:175_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [   4.37264, 92, 94, 0, 0, 1,     // 175_E_N:175_E_CA:175_E_C:176_Q_N [ 175_E_N:175_E_CA:175_E_C -- 175_E_CA:175_E_C:176_Q_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.20343, 94, 95, 0, 0, 1,     // 175_E_CA:175_E_C:176_Q_N:176_Q_CA [ 175_E_CA:175_E_C:176_Q_N -- 175_E_C:176_Q_N:176_Q_CA ] 
        [ [ 0.44758980997486436, -0.07624296316732045, 0.8909828127265595, 0.0 ], [ 0.03422519372468975, 0.9970892691065664, 0.06812947634497743, 0.0 ], [ -0.8935837946836168, 0.0, 0.4488964266719303, 15.255289178272287 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  46.46485, 95, 96, 0, 0, 1,     // 175_E_C:176_Q_N:176_Q_CA:176_Q_C [ 175_E_C:176_Q_N:176_Q_CA -- 176_Q_N:176_Q_CA:176_Q_C ] 
        [ [ -0.9907963193364446, 0.08245811146457639, 0.1073466974016873, 11.885236283375603 ], [ -0.08331467979614059, -0.9965171010652176, -0.0035116114593420397, 0.9088109362567597 ], [ 0.10668325885451983, -0.012422847430050548, 0.9942154470445072, 21.243328477742203 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 175, ' ') E sidechain

     ],
     [ // 24 : (' ', 176, ' ') Q backbone
      [ -131.31559, 96, 97, 0, 0, 1,     // 176_Q_N:176_Q_CA:176_Q_C:176_Q_O [ 176_Q_N:176_Q_CA:176_Q_C -- 176_Q_CA:176_Q_C:176_Q_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  47.59676, 96, 98, 0, 0, 1,     // 176_Q_N:176_Q_CA:176_Q_C:177_A_N [ 176_Q_N:176_Q_CA:176_Q_C -- 176_Q_CA:176_Q_C:177_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.79643, 98, 99, 0, 0, 1,     // 176_Q_CA:176_Q_C:177_A_N:177_A_CA [ 176_Q_CA:176_Q_C:177_A_N -- 176_Q_C:177_A_N:177_A_CA ] 
        [ [ 0.2903992632002363, -0.7384172131291632, 0.608611770579012, 0.0 ], [ 0.3179916630986096, 0.6743441401499385, 0.6664392566807724, 0.0 ], [ -0.9025239997543225, 0.0, 0.43063955910652213, 15.228627936918233 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.21358, 99, 100, 0, 0, 1,     // 176_Q_C:177_A_N:177_A_CA:177_A_C [ 176_Q_C:177_A_N:177_A_CA -- 177_A_N:177_A_CA:177_A_C ] 
        [ [ -0.668387932904658, 0.739444325616961, 0.08049633817954797, 8.08437228644345 ], [ -0.7347006209756016, -0.6732100763518235, 0.0836850681808731, 8.852512090231675 ], [ 0.11607139477711081, -0.0032066199103767943, 0.9932356965510485, 20.948942077038893 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 176, ' ') Q sidechain

     ],
     [ // 25 : (' ', 177, ' ') A backbone
      [ -38.45271, 100, 101, 0, 0, 1,     // 177_A_N:177_A_CA:177_A_C:177_A_O [ 177_A_N:177_A_CA:177_A_C -- 177_A_CA:177_A_C:177_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 139.92237, 100, 102, 0, 0, 1,     // 177_A_N:177_A_CA:177_A_C:178_S_N [ 177_A_N:177_A_CA:177_A_C -- 177_A_CA:177_A_C:178_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.77348, 102, 103, 0, 0, 1,     // 177_A_CA:177_A_C:178_S_N:178_S_CA [ 177_A_CA:177_A_C:178_S_N -- 177_A_C:178_S_N:178_S_CA ] 
        [ [ -0.33391198962006224, -0.6438249378584197, -0.6884709380791418, 0.0 ], [ 0.28095726801059284, -0.765172823218129, 0.5792871172056409, 0.0 ], [ -0.8997587436307554, 0.0, 0.43638767542175694, 15.228394009543013 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -75.50371, 103, 104, 0, 0, 1,     // 177_A_C:178_S_N:178_S_CA:178_S_C [ 177_A_C:178_S_N:178_S_CA -- 178_S_N:178_S_CA:178_S_C ] 
        [ [ 0.7531938994773463, 0.6508248624249082, -0.09552982905725657, -9.136846134437716 ], [ -0.6486364578236297, 0.7589835777628683, 0.05669809757205952, 7.687844126488747 ], [ 0.10940610299304146, 0.019259468732949186, 0.9938105340013259, 21.01978885972796 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 177, ' ') A sidechain

     ],
     [ // 26 : (' ', 178, ' ') S backbone
      [  -3.46788, 104, 105, 0, 0, 1,     // 178_S_N:178_S_CA:178_S_C:178_S_O [ 178_S_N:178_S_CA:178_S_C -- 178_S_CA:178_S_C:178_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.93153, 104, 106, 0, 0, 1,     // 178_S_N:178_S_CA:178_S_C:179_Q_N [ 178_S_N:178_S_CA:178_S_C -- 178_S_CA:178_S_C:179_Q_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.49542, 106, 107, 0, 0, 1,     // 178_S_CA:178_S_C:179_Q_N:179_Q_CA [ 178_S_CA:178_S_C:179_Q_N -- 178_S_C:179_Q_N:179_Q_CA ] 
        [ [ -0.4528768261591958, -0.036093855392846624, -0.8908421936184062, 0.0 ], [ 0.01635672864472669, -0.9993484045131009, 0.03217489433024607, 0.0 ], [ -0.8914230408487411, 0.0, 0.45317211106155203, 15.256456851623286 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.44947, 107, 108, 0, 0, 1,     // 178_S_C:179_Q_N:179_Q_CA:179_Q_C [ 178_S_C:179_Q_N:179_Q_CA -- 179_Q_N:179_Q_CA:179_Q_C ] 
        [ [ 0.9967023274389645, 0.04008071803017866, -0.07055498933411036, -11.834110146737697 ], [ -0.040534208252266286, 0.9991656061370817, -0.005006942588052704, 0.42741716343408626 ], [ 0.07029543682994668, 0.007850321961770658, 0.9974953253053287, 21.27647781173203 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 178, ' ') S sidechain

     ],
     [ // 27 : (' ', 179, ' ') Q backbone
      [ 133.13086, 108, 109, 0, 0, 1,     // 179_Q_N:179_Q_CA:179_Q_C:179_Q_O [ 179_Q_N:179_Q_CA:179_Q_C -- 179_Q_CA:179_Q_C:179_Q_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -46.34169, 108, 110, 0, 0, 1,     // 179_Q_N:179_Q_CA:179_Q_C:180_E_N [ 179_Q_N:179_Q_CA:179_Q_C -- 179_Q_CA:179_Q_C:180_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.51734, 110, 111, 0, 0, 1,     // 179_Q_CA:179_Q_C:180_E_N:180_E_CA [ 179_Q_CA:179_Q_C:180_E_N -- 179_Q_C:180_E_N:180_E_CA ] 
        [ [ 0.30540126769805526, 0.7234696071522584, 0.619129867810767, 0.0 ], [ -0.3200500333244761, 0.6903562323373035, -0.6488268248477747, 0.0 ], [ -0.8968266509518005, 0.0, 0.4423821403973916, 15.284172599349398 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -59.07605, 111, 112, 0, 0, 1,     // 179_Q_C:180_E_N:180_E_CA:180_E_C [ 179_Q_C:180_E_N:180_E_CA -- 180_E_N:180_E_CA:180_E_C ] 
        [ [ -0.6837261126830209, -0.7260166282826647, 0.0736101779133574, 8.241944455042823 ], [ 0.7229777196236825, -0.6876356446705699, -0.06678650392279628, -8.637274551534984 ], [ 0.09910509453658238, 0.007554841862069866, 0.9950482925975674, 21.17322628359624 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 179, ' ') Q sidechain

     ],
     [ // 28 : (' ', 180, ' ') E backbone
      [ 137.68712, 112, 113, 0, 0, 1,     // 180_E_N:180_E_CA:180_E_C:180_E_O [ 180_E_N:180_E_CA:180_E_C -- 180_E_CA:180_E_C:180_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.44142, 112, 114, 0, 0, 1,     // 180_E_N:180_E_CA:180_E_C:181_V_N [ 180_E_N:180_E_CA:180_E_C -- 180_E_CA:180_E_C:181_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.30496, 114, 115, 0, 0, 1,     // 180_E_CA:180_E_C:181_V_N:181_V_CA [ 180_E_CA:180_E_C:181_V_N -- 180_E_C:181_V_N:181_V_CA ] 
        [ [ 0.3336725153523985, 0.6748361000051233, 0.6582240428822687, 0.0 ], [ -0.3051275511861055, 0.7379676402999492, -0.6019143953896562, 0.0 ], [ -0.8919416068362177, 0.0, 0.4521506054341032, 15.221202894948204 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.96470, 115, 116, 0, 0, 1,     // 180_E_C:181_V_N:181_V_CA:181_V_C [ 180_E_C:181_V_N:181_V_CA -- 181_V_N:181_V_CA:181_V_C ] 
        [ [ -0.7315839282816977, -0.678834043106159, 0.06300236344824794, 8.760218889221694 ], [ 0.6775215754423843, -0.7342120110542653, -0.043557291395595985, -8.010801053540433 ], [ 0.07582526419333539, 0.010819646195549943, 0.997062417587893, 21.238816984449517 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 180, ' ') E sidechain

     ],
     [ // 29 : (' ', 181, ' ') V backbone
      [ 138.65051, 116, 117, 0, 0, 1,     // 181_V_N:181_V_CA:181_V_C:181_V_O [ 181_V_N:181_V_CA:181_V_C -- 181_V_CA:181_V_C:181_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -39.99620, 116, 118, 0, 0, 1,     // 181_V_N:181_V_CA:181_V_C:182_K_N [ 181_V_N:181_V_CA:181_V_C -- 181_V_CA:181_V_C:182_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.64559, 118, 119, 0, 0, 1,     // 181_V_CA:181_V_C:182_K_N:182_K_CA [ 181_V_CA:181_V_C:182_K_N -- 181_V_C:182_K_N:182_K_CA ] 
        [ [ 0.3422611220817856, 0.6427367830295172, 0.6853807351043538, 0.0 ], [ -0.2871524866526414, 0.7660870888749316, -0.5750252356795632, 0.0 ], [ -0.8946512022685275, 0.0, 0.44676529216074784, 15.268211301663726 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.72625, 119, 120, 0, 0, 1,     // 181_V_C:182_K_N:182_K_CA:182_K_C [ 181_V_C:182_K_N:182_K_CA -- 182_K_N:182_K_CA:182_K_C ] 
        [ [ -0.755802454471755, -0.6506471645654984, 0.07362687727554362, 9.106277889275185 ], [ 0.6500305537868196, -0.759085690180909, -0.035343940161982244, -7.640044899491947 ], [ 0.0788855434035375, 0.0211466830838868, 0.996659364495452, 21.204136683458344 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 181, ' ') V sidechain

     ],
     [ // 30 : (' ', 182, ' ') K backbone
      [ 137.08692, 120, 121, 0, 0, 1,     // 182_K_N:182_K_CA:182_K_C:182_K_O [ 182_K_N:182_K_CA:182_K_C -- 182_K_CA:182_K_C:182_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.94998, 120, 122, 0, 0, 1,     // 182_K_N:182_K_CA:182_K_C:183_N_N [ 182_K_N:182_K_CA:182_K_C -- 182_K_CA:182_K_C:183_N_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.12349, 122, 123, 0, 0, 1,     // 182_K_CA:182_K_C:183_N_N:183_N_CA [ 182_K_CA:182_K_C:183_N_N -- 182_K_C:183_N_N:183_N_CA ] 
        [ [ 0.3273031400529061, 0.6813595746320743, 0.6546921296065095, 0.0 ], [ -0.3046812984052789, 0.7319488575421091, -0.6094425127482285, 0.0 ], [ -0.8944506475561305, 0.0, 0.4471666793114385, 15.256263829849193 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.99990, 123, 124, 0, 0, 1,     // 182_K_C:183_N_N:183_N_CA:183_N_C [ 182_K_C:183_N_N:183_N_CA -- 183_N_N:183_N_CA:183_N_C ] 
        [ [ -0.7232489086795999, -0.6862867263648351, 0.07695157769146539, 8.733180121811356 ], [ 0.6843044674861185, -0.7272023870002525, -0.05388955482906476, -8.129578769976211 ], [ 0.09294305714956616, 0.013682746675024343, 0.9955774056150127, 21.22118532240813 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 182, ' ') K sidechain

     ],
     [ // 31 : (' ', 183, ' ') N backbone
      [ 133.63160, 124, 125, 0, 0, 1,     // 183_N_N:183_N_CA:183_N_C:183_N_O [ 183_N_N:183_N_CA:183_N_C -- 183_N_CA:183_N_C:183_N_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.61274, 124, 126, 0, 0, 1,     // 183_N_N:183_N_CA:183_N_C:184_W_N [ 183_N_N:183_N_CA:183_N_C -- 183_N_CA:183_N_C:184_W_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.17932, 126, 127, 0, 0, 1,     // 183_N_CA:183_N_C:184_W_N:184_W_CA [ 183_N_CA:183_N_C:184_W_N -- 183_N_C:184_W_N:184_W_CA ] 
        [ [ 0.31385165700007556, 0.7146282790795196, 0.6251428317897848, 0.0 ], [ -0.3206373907074687, 0.6995044122375814, -0.6386589394508323, 0.0 ], [ -0.89369390793415, 0.0, 0.44867716558945475, 15.184931855915602 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.83417, 127, 128, 0, 0, 1,     // 183_N_C:184_W_N:184_W_CA:184_W_C [ 183_N_C:184_W_N:184_W_CA -- 184_W_N:184_W_CA:184_W_C ] 
        [ [ -0.6912074239233302, -0.719050306942283, 0.07210376688829899, 8.324997854026268 ], [ 0.7169593867120545, -0.6948401281192873, -0.056271077478996165, -8.50499122762196 ], [ 0.09056232615579408, 0.012800485982216718, 0.9958085220763431, 21.15994468616563 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 183, ' ') N sidechain

     ],
     [ // 32 : (' ', 184, ' ') W backbone
      [ 135.79091, 128, 129, 0, 0, 1,     // 184_W_N:184_W_CA:184_W_C:184_W_O [ 184_W_N:184_W_CA:184_W_C -- 184_W_CA:184_W_C:184_W_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  23.53516, 128, 130, 0, 0, 1,     // 184_W_N:184_W_CA:184_W_C:186_T_N [ 184_W_N:184_W_CA:184_W_C -- 184_W_CA:184_W_C:186_T_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  74.00301, 130, 131, 0, 0, 1,     // 184_W_CA:184_W_C:186_T_N:186_T_CA [ 184_W_CA:184_W_C:186_T_N -- 184_W_C:186_T_N:186_T_CA ] 
        [ [ 0.5571352182588668, -0.3993117789821032, 0.7281143122764397, 0.0 ], [ 0.24265594218633654, 0.9168151957543832, 0.31712456630717417, 0.0 ], [ -0.7941778404723379, 0.0, 0.6076854101446684, 15.243406129108743 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  28.20627, 131, 132, 0, 0, 1,     // 184_W_C:186_T_N:186_T_CA:186_T_C [ 184_W_C:186_T_N:186_T_CA -- 186_T_N:186_T_CA:186_T_C ] 
        [ [ -0.6501813357046432, -0.6456059040745787, 0.40057115138932803, 23.796079001631245 ], [ 0.5103990131613109, 0.019402901981610793, 0.859718776553503, 10.36419845890694 ], [ -0.5628117607737831, 0.7634242228399213, 0.31690121160383095, 35.103653049308065 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 184, ' ') W sidechain

     ],
     [ // 33 : (' ', 186, ' ') T backbone
      [ 140.92472, 132, 133, 0, 0, 1,     // 186_T_N:186_T_CA:186_T_C:186_T_O [ 186_T_N:186_T_CA:186_T_C -- 186_T_CA:186_T_C:186_T_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -39.01887, 132, 134, 0, 0, 1,     // 186_T_N:186_T_CA:186_T_C:187_E_N [ 186_T_N:186_T_CA:186_T_C -- 186_T_CA:186_T_C:187_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.48606, 134, 135, 0, 0, 1,     // 186_T_CA:186_T_C:187_E_N:187_E_CA [ 186_T_CA:186_T_C:187_E_N -- 186_T_C:187_E_N:187_E_CA ] 
        [ [ 0.34533375178525055, 0.6295762715796821, 0.6959729291730746, 0.0 ], [ -0.27983409817877963, 0.7769386837220982, -0.5639673387995467, 0.0 ], [ -0.8957887459520756, 0.0, 0.44448005874910507, 15.409636813519732 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -72.50210, 135, 136, 0, 0, 1,     // 186_T_C:187_E_N:187_E_CA:187_E_C [ 186_T_C:187_E_N:187_E_CA -- 187_E_N:187_E_CA:187_E_C ] 
        [ [ -0.7680413953280812, -0.6326485474160456, 0.09933896775617204, 9.305769570950225 ], [ 0.6287193450223492, -0.7743973477969921, -0.07085713033741119, -7.540738842595829 ], [ 0.12175549374530313, 0.008035121489567757, 0.9925275999011238, 21.352725735837563 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 186, ' ') T sidechain

     ],
     [ // 34 : (' ', 187, ' ') E backbone
      [ 149.79789, 136, 137, 0, 0, 1,     // 187_E_N:187_E_CA:187_E_C:187_E_O [ 187_E_N:187_E_CA:187_E_C -- 187_E_CA:187_E_C:187_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -33.62689, 136, 138, 0, 0, 1,     // 187_E_N:187_E_CA:187_E_C:188_T_N [ 187_E_N:187_E_CA:187_E_C -- 187_E_CA:187_E_C:188_T_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -177.60466, 138, 139, 0, 0, 1,     // 187_E_CA:187_E_C:188_T_N:188_T_CA [ 187_E_CA:187_E_C:188_T_N -- 187_E_C:188_T_N:188_T_CA ] 
        [ [ 0.3704938651428755, 0.5537823297307838, 0.7456939232482978, 0.0 ], [ -0.24640620714330716, 0.832661474596937, -0.4959423855684288, 0.0 ], [ -0.8955547314222299, 0.0, 0.44495137153093217, 15.204033491582353 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -100.08243, 139, 140, 0, 0, 1,     // 187_E_C:188_T_N:188_T_CA:188_T_C [ 187_E_C:188_T_N:188_T_CA -- 188_T_N:188_T_CA:188_T_C ] 
        [ [ -0.8408380380774797, -0.5378138817978889, 0.06121782638661146, 9.927045710604874 ], [ 0.5326414464199848, -0.8422323225996899, -0.08329348308216197, -6.602229920177184 ], [ 0.09635602356698467, -0.03742917730619258, 0.9946429376457364, 21.127445854851842 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 187, ' ') E sidechain

     ],
     [ // 35 : (' ', 188, ' ') T backbone
      [  90.06146, 140, 141, 0, 0, 1,     // 188_T_N:188_T_CA:188_T_C:188_T_O [ 188_T_N:188_T_CA:188_T_C -- 188_T_CA:188_T_C:188_T_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -89.69692, 140, 142, 0, 0, 1,     // 188_T_N:188_T_CA:188_T_C:189_L_N [ 188_T_N:188_T_CA:188_T_C -- 188_T_CA:188_T_C:189_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -177.75783, 142, 143, 0, 0, 1,     // 188_T_CA:188_T_C:189_L_N:189_L_CA [ 188_T_CA:188_T_C:189_L_N -- 188_T_C:189_L_N:189_L_CA ] 
        [ [ 0.0023304302256365976, 0.9999860095636259, 0.004748660019338789, 0.0 ], [ -0.4405557959938835, 0.005289676456630151, -0.8977096467896398, 0.0 ], [ -0.8977222063150488, 0.0, 0.4405619596479488, 15.232859526498055 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -53.41239, 143, 144, 0, 0, 1,     // 188_T_C:189_L_N:189_L_CA:189_L_C [ 188_T_C:189_L_N:189_L_CA -- 189_L_N:189_L_CA:189_L_C ] 
        [ [ -0.02629541845040355, -0.9991292386123587, -0.032393139989918905, 0.06319527295238384 ], [ 0.9934745042660834, -0.022521603686339674, -0.11180870601458161, -11.946739907642131 ], [ 0.1109818018495729, -0.035121815404158103, 0.993201640021269, 21.095867865483317 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 188, ' ') T sidechain

     ],
     [ // 36 : (' ', 189, ' ') L backbone
      [ 139.85449, 144, 145, 0, 0, 1,     // 189_L_N:189_L_CA:189_L_C:189_L_O [ 189_L_N:189_L_CA:189_L_C -- 189_L_CA:189_L_C:189_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -39.58667, 144, 146, 0, 0, 1,     // 189_L_N:189_L_CA:189_L_C:190_L_N [ 189_L_N:189_L_CA:189_L_C -- 189_L_CA:189_L_C:190_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.97719, 146, 147, 0, 0, 1,     // 189_L_CA:189_L_C:190_L_N:190_L_CA [ 189_L_CA:189_L_C:190_L_N -- 189_L_C:190_L_N:190_L_CA ] 
        [ [ 0.36561188287973495, 0.6372446712506888, 0.6784151973973876, 0.0 ], [ -0.3023171757626279, 0.7706615527990232, -0.5609679993335688, 0.0 ], [ -0.8803023777862032, 0.0, 0.4744130306641636, 15.238268545680574 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.30473, 147, 148, 0, 0, 1,     // 189_L_C:190_L_N:190_L_CA:190_L_C [ 189_L_C:190_L_N:190_L_CA -- 190_L_N:190_L_CA:190_L_C ] 
        [ [ -0.7703612060235341, -0.6373901858105528, 0.01665422728182822, 9.03452081829112 ], [ 0.6372510548690136, -0.7705411268667296, -0.013321594322001187, -7.470465119025887 ], [ 0.02132382053702769, 0.00035048443532023854, 0.9997725600547183, 21.556072773013284 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 189, ' ') L sidechain

     ],
     [ // 37 : (' ', 190, ' ') L backbone
      [ 140.35392, 148, 149, 0, 0, 1,     // 190_L_N:190_L_CA:190_L_C:190_L_O [ 190_L_N:190_L_CA:190_L_C -- 190_L_CA:190_L_C:190_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -39.28926, 148, 150, 0, 0, 1,     // 190_L_N:190_L_CA:190_L_C:191_V_N [ 190_L_N:190_L_CA:190_L_C -- 190_L_CA:190_L_C:191_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.39889, 150, 151, 0, 0, 1,     // 190_L_CA:190_L_C:191_V_N:191_V_CA [ 190_L_CA:190_L_C:191_V_N -- 190_L_C:191_V_N:191_V_CA ] 
        [ [ 0.3396544097582195, 0.633235854800424, 0.6954475063777078, 0.0 ], [ -0.2778976445549742, 0.7739588827548763, -0.5689996018940794, 0.0 ], [ -0.8985587243372537, 0.0, 0.43885329999603195, 15.253068570552784 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -75.70206, 151, 152, 0, 0, 1,     // 190_L_C:191_V_N:191_V_CA:191_V_C [ 190_L_C:191_V_N:191_V_CA -- 191_V_N:191_V_CA:191_V_C ] 
        [ [ -0.7659862982354297, -0.6367643998600637, 0.08829546979565957, 9.244423407891194 ], [ 0.633928424775217, -0.7710007982375354, -0.0607661203228852, -7.563580558693635 ], [ 0.10676957983243814, 0.009426992518094672, 0.9942391003347577, 21.086644335666286 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 190, ' ') L sidechain

     ],
     [ // 38 : (' ', 191, ' ') V backbone
      [ 130.04069, 152, 153, 0, 0, 1,     // 191_V_N:191_V_CA:191_V_C:191_V_O [ 191_V_N:191_V_CA:191_V_C -- 191_V_CA:191_V_C:191_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -47.86944, 152, 154, 0, 0, 1,     // 191_V_N:191_V_CA:191_V_C:192_Q_N [ 191_V_N:191_V_CA:191_V_C -- 191_V_CA:191_V_C:192_Q_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.07306, 154, 155, 0, 0, 1,     // 191_V_CA:191_V_C:192_Q_N:192_Q_CA [ 191_V_CA:191_V_C:192_Q_N -- 191_V_C:192_Q_N:192_Q_CA ] 
        [ [ 0.2946050315448642, 0.7416181093688072, 0.602669441107383, 0.0 ], [ -0.32569641951688094, 0.6708223161584861, -0.6662726637478443, 0.0 ], [ -0.8984039835743903, 0.0, 0.4391699924831688, 15.310336220149718 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.71576, 155, 156, 0, 0, 1,     // 191_V_C:192_Q_N:192_Q_CA:192_Q_C [ 191_V_C:192_Q_N:192_Q_CA -- 192_Q_N:192_Q_CA:192_Q_C ] 
        [ [ -0.672706381415069, -0.7367550901537734, 0.06824999293735434, 8.018109356917229 ], [ 0.7306533123457872, -0.6760034781173178, -0.09573418789240636, -8.864307222278558 ], [ 0.11666988283860204, -0.014533915707556454, 0.9930644005967753, 21.15319598011412 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 191, ' ') V sidechain

     ],
     [ // 39 : (' ', 192, ' ') Q backbone
      [ 154.70540, 156, 157, 0, 0, 1,     // 192_Q_N:192_Q_CA:192_Q_C:192_Q_O [ 192_Q_N:192_Q_CA:192_Q_C -- 192_Q_CA:192_Q_C:192_Q_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -26.31002, 156, 158, 0, 0, 1,     // 192_Q_N:192_Q_CA:192_Q_C:193_N_N [ 192_Q_N:192_Q_CA:192_Q_C -- 192_Q_CA:192_Q_C:193_N_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.99302, 158, 159, 0, 0, 1,     // 192_Q_CA:192_Q_C:193_N_N:193_N_CA [ 192_Q_CA:192_Q_C:193_N_N -- 192_Q_C:193_N_N:193_N_CA ] 
        [ [ 0.40769205607459125, 0.4432279880766864, 0.7983333501734504, 0.0 ], [ -0.20158269945773694, 0.8964089192915768, -0.3947346763255528, 0.0 ], [ -0.8905905920752835, 0.0, 0.45480588970130537, 15.258649771447994 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -97.39357, 159, 160, 0, 0, 1,     // 192_Q_C:193_N_N:193_N_CA:193_N_C [ 192_Q_C:193_N_N:193_N_CA -- 193_N_N:193_N_CA:193_N_C ] 
        [ [ -0.8940215313295972, -0.4432776672402091, 0.06504161164331916, 10.612379173125152 ], [ 0.4421179923592498, -0.8963843472235371, -0.03204345307328248, -5.247274283402357 ], [ 0.07250642972390922, 0.00010852977389018585, 0.9973679390626009, 21.30446077337111 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 192, ' ') Q sidechain

     ],
     [ // 40 : (' ', 193, ' ') N backbone
      [ -174.25901, 160, 161, 0, 0, 1,     // 193_N_N:193_N_CA:193_N_C:193_N_O [ 193_N_N:193_N_CA:193_N_C -- 193_N_CA:193_N_C:193_N_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [   5.01268, 160, 162, 0, 0, 1,     // 193_N_N:193_N_CA:193_N_C:194_A_N [ 193_N_N:193_N_CA:193_N_C -- 193_N_CA:193_N_C:194_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.10894, 162, 163, 0, 0, 1,     // 193_N_CA:193_N_C:194_A_N:194_A_CA [ 193_N_CA:193_N_C:194_A_N -- 193_N_C:194_A_N:194_A_CA ] 
        [ [ 0.4449591611389792, -0.08737618007953937, 0.8912781541545849, 0.0 ], [ 0.039028099137858485, 0.9961753877479146, 0.07817547136396162, 0.0 ], [ -0.8947000348698895, 0.0, 0.44666749109804105, 15.266587400207582 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -76.36918, 163, 164, 0, 0, 1,     // 193_N_C:194_A_N:194_A_CA:194_A_C [ 193_N_C:194_A_N:194_A_CA -- 194_A_N:194_A_CA:194_A_C ] 
        [ [ -0.9920803576611398, 0.08044591339643332, 0.09646252619934133, 11.861178300958537 ], [ -0.07870961362413358, -0.9966618595223558, 0.021677972612304645, 1.0403634379314122 ], [ 0.09788442504343918, 0.013913762656225565, 0.995100520823229, 21.210862740582165 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 193, ' ') N sidechain

     ],
     [ // 41 : (' ', 194, ' ') A backbone
      [ -21.92381, 164, 165, 0, 0, 1,     // 194_A_N:194_A_CA:194_A_C:194_A_O [ 194_A_N:194_A_CA:194_A_C -- 194_A_CA:194_A_C:194_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 157.52778, 164, 166, 0, 0, 1,     // 194_A_N:194_A_CA:194_A_C:195_N_N [ 194_A_N:194_A_CA:194_A_C -- 194_A_CA:194_A_C:195_N_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.29573, 166, 167, 0, 0, 1,     // 194_A_CA:194_A_C:195_N_N:195_N_CA [ 194_A_CA:194_A_C:195_N_N -- 194_A_C:195_N_N:195_N_CA ] 
        [ [ -0.4233489364168933, -0.3822354612933217, -0.8213840333026731, 0.0 ], [ 0.17511645049926164, -0.9240649609902334, 0.3397619411215124, 0.0 ], [ -0.8888812669863309, 0.0, 0.45813763565196775, 15.35185500388245 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -69.24681, 167, 168, 0, 0, 1,     // 194_A_C:195_N_N:195_N_CA:195_N_C [ 194_A_C:195_N_N:195_N_CA -- 195_N_N:195_N_CA:195_N_C ] 
        [ [ 0.9250292259155171, 0.3770030105934443, -0.046793815890764576, -10.942519133025591 ], [ -0.3758994101984588, 0.9261475909657573, 0.03082650094927006, 4.526325556204314 ], [ 0.054959663523263265, -0.010925646516513565, 0.9984287984796942, 21.455187477394812 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 194, ' ') A sidechain

     ],
     [ // 42 : (' ', 195, ' ') N backbone
      [ -12.45715, 168, 169, 0, 0, 1,     // 195_N_N:195_N_CA:195_N_C:195_N_O [ 195_N_N:195_N_CA:195_N_C -- 195_N_CA:195_N_C:195_N_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 165.14971, 168, 170, 0, 0, 1,     // 195_N_N:195_N_CA:195_N_C:196_P_N [ 195_N_N:195_N_CA:195_N_C -- 195_N_CA:195_N_C:196_P_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.73907, 170, 171, 0, 0, 1,     // 195_N_CA:195_N_C:196_P_N:196_P_CA [ 195_N_CA:195_N_C:196_P_N -- 195_N_C:196_P_N:196_P_CA ] 
        [ [ -0.4502434344010165, -0.2562943494098208, -0.8553327166894612, 0.0 ], [ 0.11938236440704063, -0.9665987825673052, 0.22679207351225456, 0.0 ], [ -0.8848890895741466, 0.0, 0.46580178096765357, 15.342454720948151 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.16922, 171, 172, 0, 0, 1,     // 195_N_C:196_P_N:196_P_CA:196_P_C [ 195_N_C:196_P_N:196_P_CA -- 196_P_N:196_P_CA:196_P_C ] 
        [ [ 0.9625721836093926, 0.2542412146903941, -0.09389460098535436, -11.510947999949428 ], [ -0.2526334768565516, 0.9671324441380802, 0.028829878023910002, 3.0521359864550544 ], [ 0.09813825815051033, -0.00402991913767882, 0.9951646306209478, 21.611149901643742 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 195, ' ') N sidechain

     ],
     [ // 43 : (' ', 196, ' ') P backbone
      [ 128.72143, 172, 173, 0, 0, 1,     // 196_P_N:196_P_CA:196_P_C:196_P_O [ 196_P_N:196_P_CA:196_P_C -- 196_P_CA:196_P_C:196_P_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -50.91642, 172, 174, 0, 0, 1,     // 196_P_N:196_P_CA:196_P_C:197_D_N [ 196_P_N:196_P_CA:196_P_C -- 196_P_CA:196_P_C:197_D_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.90236, 174, 175, 0, 0, 1,     // 196_P_CA:196_P_C:197_D_N:197_D_CA [ 196_P_CA:196_P_C:197_D_N -- 196_P_C:197_D_N:197_D_CA ] 
        [ [ 0.28328144144150086, 0.7762271706109342, 0.5632255361222261, 0.0 ], [ -0.34878197531285016, 0.630453312787984, -0.6934549402026945, 0.0 ], [ -0.893365971274758, 0.0, 0.44932976906533595, 15.297546571096374 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -59.40454, 175, 176, 0, 0, 1,     // 196_P_C:197_D_N:197_D_CA:197_D_C [ 196_P_C:197_D_N:197_D_CA -- 197_D_N:197_D_CA:197_D_C ] 
        [ [ -0.6286539859574979, -0.7757432735458188, 0.05492303240141481, 7.483988366708356 ], [ 0.7725873254799249, -0.6310467935144244, -0.06991973185689254, -9.214441413727304 ], [ 0.08889876516311335, -0.001522479418658623, 0.9960395030363471, 21.268118835367726 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 196, ' ') P sidechain

     ],
     [ // 44 : (' ', 197, ' ') D backbone
      [ 120.81353, 176, 177, 0, 0, 1,     // 197_D_N:197_D_CA:197_D_C:197_D_O [ 197_D_N:197_D_CA:197_D_C -- 197_D_CA:197_D_C:197_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -58.02023, 176, 178, 0, 0, 1,     // 197_D_N:197_D_CA:197_D_C:198_C_N [ 197_D_N:197_D_CA:197_D_C -- 197_D_CA:197_D_C:198_C_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -178.99418, 178, 179, 0, 0, 1,     // 197_D_CA:197_D_C:198_C_N:198_C_CA [ 197_D_CA:197_D_C:198_C_N -- 197_D_C:198_C_N:198_C_CA ] 
        [ [ 0.2323695567581231, 0.8482351377243935, 0.4759217795202406, 0.0 ], [ -0.37216134376485377, 0.5296198175386558, -0.7622327617450323, 0.0 ], [ -0.8986102176690246, 0.0, 0.4387478509358516, 15.231209373274181 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -67.54270, 179, 180, 0, 0, 1,     // 197_D_C:198_C_N:198_C_CA:198_C_C [ 197_D_C:198_C_N:198_C_CA -- 198_C_N:198_C_CA:198_C_C ] 
        [ [ -0.5346957872194145, -0.8440254377317469, 0.04149066872856751, 6.322943553478974 ], [ 0.8388562410862876, -0.5360711114899855, -0.09459371129454469, -10.126779093792775 ], [ 0.10208144748383453, -0.015774152516840063, 0.994650971040586, 21.060272078650144 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 197, ' ') D sidechain

     ],
     [ // 45 : (' ', 198, ' ') C backbone
      [ 142.43101, 180, 181, 0, 0, 1,     // 198_C_N:198_C_CA:198_C_C:198_C_O [ 198_C_N:198_C_CA:198_C_C -- 198_C_CA:198_C_C:198_C_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -35.81689, 180, 182, 0, 0, 1,     // 198_C_N:198_C_CA:198_C_C:199_K_N [ 198_C_N:198_C_CA:198_C_C -- 198_C_CA:198_C_C:199_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -178.97019, 182, 183, 0, 0, 1,     // 198_C_CA:198_C_C:199_K_N:199_K_CA [ 198_C_CA:198_C_C:199_K_N -- 198_C_C:199_K_N:199_K_CA ] 
        [ [ 0.36846405790723663, 0.5851966780880659, 0.7223427759625103, 0.0 ], [ -0.2659097692766729, 0.810891390973349, -0.5212937239681041, 0.0 ], [ -0.8908008939340818, 0.0, 0.4543938460919567, 15.316950958956747 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -57.33163, 183, 184, 0, 0, 1,     // 198_C_C:199_K_N:199_K_CA:199_K_C [ 198_C_C:199_K_N:199_K_CA -- 199_K_N:199_K_CA:199_K_C ] 
        [ [ -0.8134705805976722, -0.5784798862528159, 0.060221555136141955, 9.609260613332715 ], [ 0.5753409980041874, -0.8155395145641918, -0.06227388055957399, -6.93472325936879 ], [ 0.08513724518465372, -0.016010040141474634, 0.9962405974949207, 21.361711193976838 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 198, ' ') C sidechain

     ],
     [ // 46 : (' ', 199, ' ') K backbone
      [ 139.26011, 184, 185, 0, 0, 1,     // 199_K_N:199_K_CA:199_K_C:199_K_O [ 199_K_N:199_K_CA:199_K_C -- 199_K_CA:199_K_C:199_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -40.32676, 184, 186, 0, 0, 1,     // 199_K_N:199_K_CA:199_K_C:200_T_N [ 199_K_N:199_K_CA:199_K_C -- 199_K_CA:199_K_C:200_T_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.67963, 186, 187, 0, 0, 1,     // 199_K_CA:199_K_C:200_T_N:200_T_CA [ 199_K_CA:199_K_C:200_T_N -- 199_K_C:200_T_N:200_T_CA ] 
        [ [ 0.34150248409470296, 0.6471459589852593, 0.6815997073988487, 0.0 ], [ -0.28988952350594105, 0.7623661244894404, -0.5785861702395326, 0.0 ], [ -0.8940582293780678, 0.0, 0.4479507589918276, 15.270056476614418 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -65.46305, 187, 188, 0, 0, 1,     // 199_K_C:200_T_N:200_T_CA:200_T_C [ 199_K_C:200_T_N:200_T_CA -- 200_T_N:200_T_CA:200_T_C ] 
        [ [ -0.7574262525356352, -0.6490453201709772, 0.07103269904753721, 9.06956796167427 ], [ 0.6468094465349085, -0.760733318569197, -0.054058837315482465, -7.69883926842945 ], [ 0.08912357624685811, 0.00499903819251883, 0.9960080309787267, 21.230622198780416 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 199, ' ') K sidechain

     ],
     [ // 47 : (' ', 200, ' ') T backbone
      [ 132.57189, 188, 189, 0, 0, 1,     // 200_T_N:200_T_CA:200_T_C:200_T_O [ 200_T_N:200_T_CA:200_T_C -- 200_T_CA:200_T_C:200_T_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -46.34314, 188, 190, 0, 0, 1,     // 200_T_N:200_T_CA:200_T_C:201_I_N [ 200_T_N:200_T_CA:200_T_C -- 200_T_CA:200_T_C:201_I_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.35264, 190, 191, 0, 0, 1,     // 200_T_CA:200_T_C:201_I_N:201_I_CA [ 200_T_CA:200_T_C:201_I_N -- 200_T_C:201_I_N:201_I_CA ] 
        [ [ 0.2961851929191512, 0.7234870913485083, 0.6235709744267451, 0.0 ], [ -0.31040764373614327, 0.6903379090359121, -0.6535140901756444, 0.0 ], [ -0.9032836908776892, 0.0, 0.42904379006621135, 15.281478227695043 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -59.24536, 191, 192, 0, 0, 1,     // 200_T_C:201_I_N:201_I_CA:201_I_C [ 200_T_C:201_I_N:201_I_CA -- 201_I_N:201_I_CA:201_I_C ] 
        [ [ -0.6799563210518306, -0.7267873270926608, 0.09715751457898458, 8.278970615704225 ], [ 0.7215019580304124, -0.6867867397551115, -0.08808461077183576, -8.676516661935231 ], [ 0.13074527150127044, 0.01020564912440583, 0.9913634644801111, 20.977767946878913 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 200, ' ') T sidechain

     ],
     [ // 48 : (' ', 201, ' ') I backbone
      [ 133.16784, 192, 193, 0, 0, 1,     // 201_I_N:201_I_CA:201_I_C:201_I_O [ 201_I_N:201_I_CA:201_I_C -- 201_I_CA:201_I_C:201_I_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.31830, 192, 194, 0, 0, 1,     // 201_I_N:201_I_CA:201_I_C:202_L_N [ 201_I_N:201_I_CA:201_I_C -- 201_I_CA:201_I_C:202_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.17707, 194, 195, 0, 0, 1,     // 201_I_CA:201_I_C:202_L_N:202_L_CA [ 201_I_CA:201_I_C:202_L_N -- 201_I_C:202_L_N:202_L_CA ] 
        [ [ 0.3110147175508644, 0.7110240495463771, 0.6306462133664382, 0.0 ], [ -0.3144896259148309, 0.7031676905025366, -0.6376923037212215, 0.0 ], [ -0.8968646055334695, 0.0, 0.44230518801082835, 15.273603130105672 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.74651, 195, 196, 0, 0, 1,     // 201_I_C:202_L_N:202_L_CA:202_L_C [ 201_I_C:202_L_N:202_L_CA -- 202_L_N:202_L_CA:202_L_C ] 
        [ [ -0.7050161592302169, -0.7064838273431102, 0.06190974823809234, 8.386293651077848 ], [ 0.7020668792845152, -0.7076119562825007, -0.079921313413596, -8.479992117119227 ], [ 0.10027119345043917, -0.012881033680680361, 0.9948767595714264, 21.15534968060578 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 201, ' ') I sidechain

     ],
     [ // 49 : (' ', 202, ' ') L backbone
      [ 143.20176, 196, 197, 0, 0, 1,     // 202_L_N:202_L_CA:202_L_C:202_L_O [ 202_L_N:202_L_CA:202_L_C -- 202_L_CA:202_L_C:202_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -38.24054, 196, 198, 0, 0, 1,     // 202_L_N:202_L_CA:202_L_C:203_K_N [ 202_L_N:202_L_CA:202_L_C -- 202_L_CA:202_L_C:203_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.20123, 198, 199, 0, 0, 1,     // 202_L_CA:202_L_C:203_K_N:203_K_CA [ 202_L_CA:202_L_C:203_K_N -- 202_L_C:203_K_N:203_K_CA ] 
        [ [ 0.3559554195450124, 0.6189642951661651, 0.7001278030516904, 0.0 ], [ -0.2805173549568576, 0.7854191246140193, -0.5517488670206364, 0.0 ], [ -0.8914066147749536, 0.0, 0.45320442091340796, 15.234215053048274 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -67.41298, 199, 200, 0, 0, 1,     // 202_L_C:203_K_N:203_K_CA:203_K_C [ 202_L_C:203_K_N:203_K_CA -- 203_K_N:203_K_CA:203_K_C ] 
        [ [ -0.778913678571602, -0.6238664122986812, 0.06390759688508908, 9.300763977655862 ], [ 0.6229987301845298, -0.7814321938407353, -0.03516118055070235, -7.329641766446315 ], [ 0.07187533319935734, 0.012426827222887474, 0.9973362073256209, 21.25475492330194 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 202, ' ') L sidechain

     ],
     [ // 50 : (' ', 203, ' ') K backbone
      [ 144.47065, 200, 201, 0, 0, 1,     // 203_K_N:203_K_CA:203_K_C:203_K_O [ 203_K_N:203_K_CA:203_K_C -- 203_K_CA:203_K_C:203_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -32.92679, 200, 202, 0, 0, 1,     // 203_K_N:203_K_CA:203_K_C:204_A_N [ 203_K_N:203_K_CA:203_K_C -- 203_K_CA:203_K_C:204_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.84042, 202, 203, 0, 0, 1,     // 203_K_CA:203_K_C:204_A_N:204_A_CA [ 203_K_CA:203_K_C:204_A_N -- 203_K_C:204_A_N:204_A_CA ] 
        [ [ 0.3633625202278875, 0.5435670392794, 0.7566389843926, 0.0 ], [ -0.235310873344971, 0.8393657568718342, -0.48999379493683365, 0.0 ], [ -0.9014413301925229, 0.0, 0.4329012915443138, 15.284294617626834 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -74.78236, 203, 204, 0, 0, 1,     // 203_K_C:204_A_N:204_A_CA:204_A_C [ 203_K_C:204_A_N:204_A_CA -- 204_A_N:204_A_CA:204_A_C ] 
        [ [ -0.8331625054985081, -0.5425529123219108, 0.10712878587224453, 10.082370316601079 ], [ 0.537739176845321, -0.8400178771764919, -0.07215638370661771, -6.52926824455883 ], [ 0.12913875141554054, -0.002510648305337907, 0.9916233556789201, 21.052793494787444 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 203, ' ') K sidechain

     ],
     [ // 51 : (' ', 204, ' ') A backbone
      [ -172.71461, 204, 205, 0, 0, 1,     // 204_A_N:204_A_CA:204_A_C:204_A_O [ 204_A_N:204_A_CA:204_A_C -- 204_A_CA:204_A_C:204_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [   5.51567, 204, 206, 0, 0, 1,     // 204_A_N:204_A_CA:204_A_C:205_L_N [ 204_A_N:204_A_CA:204_A_C -- 204_A_CA:204_A_C:205_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.23460, 206, 207, 0, 0, 1,     // 204_A_CA:204_A_C:205_L_N:205_L_CA [ 204_A_CA:204_A_C:205_L_N -- 204_A_C:205_L_N:205_L_CA ] 
        [ [ 0.4518216594334419, -0.09611793517070269, 0.8869151766687339, 0.0 ], [ 0.04363017474730206, 0.9953699525997962, 0.08564499584330243, 0.0 ], [ -0.8910407375189592, 0.0, 0.45392334604167395, 15.350555032942154 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -90.49941, 207, 208, 0, 0, 1,     // 204_A_C:205_L_N:205_L_CA:205_L_C [ 204_A_C:205_L_N:205_L_CA -- 205_L_N:205_L_CA:205_L_C ] 
        [ [ -0.993043555815724, 0.09007379031589617, 0.07583672297107188, 11.805822680003057 ], [ -0.0888958090669965, -0.9958639645708529, 0.018774962038500754, 1.1400296904979907 ], [ 0.07721419159187794, 0.011902788217514908, 0.9969434749520465, 21.392777192010307 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 204, ' ') A sidechain

     ],
     [ // 52 : (' ', 205, ' ') L backbone
      [ 153.61219, 208, 209, 0, 0, 1,     // 205_L_N:205_L_CA:205_L_C:205_L_O [ 205_L_N:205_L_CA:205_L_C -- 205_L_CA:205_L_C:205_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -25.34807, 208, 210, 0, 0, 1,     // 205_L_N:205_L_CA:205_L_C:206_G_N [ 205_L_N:205_L_CA:205_L_C -- 205_L_CA:205_L_C:206_G_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.90507, 210, 211, 0, 0, 1,     // 205_L_CA:205_L_C:206_G_N:206_G_CA [ 205_L_CA:205_L_C:206_G_N -- 205_L_C:206_G_N:206_G_CA ] 
        [ [ 0.4007671507915442, 0.4281162212230028, 0.8100013530681073, 0.0 ], [ -0.18985329342789117, 0.9037236862701657, -0.3837176381779265, 0.0 ], [ -0.8962931539518812, 0.0, 0.4434620414184164, 15.26583503186085 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  84.64211, 211, 212, 0, 0, 1,     // 205_L_C:206_G_N:206_G_CA:206_G_C [ 205_L_C:206_G_N:206_G_CA -- 206_G_N:206_G_CA:206_G_C ] 
        [ [ -0.9014085268356037, -0.42745165893756115, 0.06890389698408732, 10.782541057437834 ], [ 0.42608251326210983, -0.9040369871960857, -0.034217213119543774, -5.107955897165864 ], [ 0.07691787594773282, -0.001484942069974408, 0.9970363259714954, 21.16909389790773 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 205, ' ') L sidechain

     ],
     [ // 53 : (' ', 206, ' ') G backbone
      [  -4.33158, 212, 213, 0, 0, 1,     // 206_G_N:206_G_CA:206_G_C:206_G_O [ 206_G_N:206_G_CA:206_G_C -- 206_G_CA:206_G_C:206_G_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.27810, 212, 214, 0, 0, 1,     // 206_G_N:206_G_CA:206_G_C:207_P_N [ 206_G_N:206_G_CA:206_G_C -- 206_G_CA:206_G_C:207_P_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.95635, 214, 215, 0, 0, 1,     // 206_G_CA:206_G_C:207_P_N:207_P_CA [ 206_G_CA:206_G_C:207_P_N -- 206_G_C:207_P_N:207_P_CA ] 
        [ [ -0.4718692807282836, -0.08231949445852135, -0.8778171123514656, 0.0 ], [ 0.038976326652597086, -0.9966059907667089, 0.0725075522175213, 0.0 ], [ -0.8808065780099751, 0.0, 0.4734762635384774, 15.14728937361982 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -65.34438, 215, 216, 0, 0, 1,     // 206_G_C:207_P_N:207_P_CA:207_P_C [ 206_G_C:207_P_N:207_P_CA -- 207_P_N:207_P_CA:207_P_C ] 
        [ [ 0.9942456658270479, 0.08196001705294768, -0.06897906630863238, -11.781911133758783 ], [ -0.08171686101229815, 0.996635392411429, 0.006344227235335269, 0.9731839636445476 ], [ 0.06926695179107366, -0.0006709676574533954, 0.9975979346369838, 21.502206796009915 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 206, ' ') G sidechain

     ],
     [ // 54 : (' ', 207, ' ') P backbone
      [ -25.40292, 216, 217, 0, 0, 1,     // 207_P_N:207_P_CA:207_P_C:207_P_O [ 207_P_N:207_P_CA:207_P_C -- 207_P_CA:207_P_C:207_P_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 152.97004, 216, 218, 0, 0, 1,     // 207_P_N:207_P_CA:207_P_C:208_G_N [ 207_P_N:207_P_CA:207_P_C -- 207_P_CA:207_P_C:208_G_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.93518, 218, 219, 0, 0, 1,     // 207_P_CA:207_P_C:208_G_N:208_G_CA [ 207_P_CA:207_P_C:208_G_N -- 207_P_C:208_G_N:208_G_CA ] 
        [ [ -0.38100431828626524, -0.454456339675305, -0.8051739841650014, 0.0 ], [ 0.19438241033728948, -0.8907690134535011, 0.41078710206449015, 0.0 ], [ -0.9039088383231375, 0.0, 0.427725159420528, 15.31825729298767 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  78.00817, 219, 220, 0, 0, 1,     // 207_P_C:208_G_N:208_G_CA:208_G_C [ 207_P_C:208_G_N:208_G_CA -- 208_G_N:208_G_CA:208_G_C ] 
        [ [ 0.8818735702094105, 0.4614582674213084, -0.09672266329373, -10.70599194142724 ], [ -0.46068960747225207, 0.8870028751957435, 0.03147991362038239, 5.4620286929730675 ], [ 0.1003199468359624, 0.016797821972168805, 0.9948134204180304, 21.005502782731973 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 207, ' ') P sidechain

     ],
     [ // 55 : (' ', 208, ' ') G backbone
      [ -176.60769, 220, 221, 0, 0, 1,     // 208_G_N:208_G_CA:208_G_C:208_G_O [ 208_G_N:208_G_CA:208_G_C -- 208_G_CA:208_G_C:208_G_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [   2.35721, 220, 222, 0, 0, 1,     // 208_G_N:208_G_CA:208_G_C:209_A_N [ 208_G_N:208_G_CA:208_G_C -- 208_G_CA:208_G_C:209_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.08404, 222, 223, 0, 0, 1,     // 208_G_CA:208_G_C:209_A_N:209_A_CA [ 208_G_CA:208_G_C:209_A_N -- 208_G_C:209_A_N:209_A_CA ] 
        [ [ 0.4664881594190236, -0.04112950751217321, 0.8835706880231246, 0.0 ], [ 0.019202677105530204, 0.9991538237988213, 0.03637160403631738, 0.0 ], [ -0.8843189776963019, 0.0, 0.4668832248926573, 15.161626375823362 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -65.16037, 223, 224, 0, 0, 1,     // 208_G_C:209_A_N:209_A_CA:209_A_C [ 208_G_C:209_A_N:209_A_CA -- 209_A_N:209_A_CA:209_A_C ] 
        [ [ -0.9978210461728408, 0.033667050811137786, 0.056742307885896036, 11.767913716999987 ], [ -0.03280633209819542, -0.9993331211166279, 0.01603301635881408, 0.4844183989465587 ], [ 0.05724425201538113, 0.014136574159935376, 0.9982601128374409, 21.379851562473846 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 208, ' ') G sidechain

     ],
     [ // 56 : (' ', 209, ' ') A backbone
      [ -47.87201, 224, 225, 0, 0, 1,     // 209_A_N:209_A_CA:209_A_C:209_A_O [ 209_A_N:209_A_CA:209_A_C -- 209_A_CA:209_A_C:209_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 130.56362, 224, 226, 0, 0, 1,     // 209_A_N:209_A_CA:209_A_C:210_T_N [ 209_A_N:209_A_CA:209_A_C -- 209_A_CA:209_A_C:210_T_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.18064, 226, 227, 0, 0, 1,     // 209_A_CA:209_A_C:210_T_N:210_T_CA [ 209_A_CA:209_A_C:210_T_N -- 209_A_C:210_T_N:210_T_CA ] 
        [ [ -0.29036814395454785, -0.7596843366202122, -0.5818642880176605, 0.0 ], [ 0.3392139603514159, -0.6502920180149131, 0.6797456733284882, 0.0 ], [ -0.8947738429788272, 0.0, 0.44651961874132806, 15.266424325693473 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -78.96673, 227, 228, 0, 0, 1,     // 209_A_C:210_T_N:210_T_CA:210_T_C [ 209_A_C:210_T_N:210_T_CA -- 210_T_N:210_T_CA:210_T_C ] 
        [ [ 0.6538239846648356, 0.755454382065359, -0.042460260188128554, -7.734842626879575 ], [ -0.752457051510823, 0.6550762985027251, 0.06843558118119031, 9.036000176966642 ], [ 0.0795146697700217, -0.012795302193198079, 0.9967515726263737, 21.20210220093302 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 209, ' ') A sidechain

     ],
     [ // 57 : (' ', 210, ' ') T backbone
      [ -17.07508, 228, 229, 0, 0, 1,     // 210_T_N:210_T_CA:210_T_C:210_T_O [ 210_T_N:210_T_CA:210_T_C -- 210_T_CA:210_T_C:210_T_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 162.73827, 228, 230, 0, 0, 1,     // 210_T_N:210_T_CA:210_T_C:211_L_N [ 210_T_N:210_T_CA:210_T_C -- 210_T_CA:210_T_C:211_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.31591, 230, 231, 0, 0, 1,     // 210_T_CA:210_T_C:211_L_N:211_L_CA [ 210_T_CA:210_T_C:211_L_N -- 210_T_C:211_L_N:211_L_CA ] 
        [ [ -0.41489615484416426, -0.29673703805498997, -0.8601211024860852, 0.0 ], [ 0.12892179296196965, -0.9549592296251979, 0.26726773290676076, 0.0 ], [ -0.9006888208449122, 0.0, 0.43446478338871347, 15.21184516945889 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -57.48376, 231, 232, 0, 0, 1,     // 210_T_C:211_L_N:211_L_CA:211_L_C [ 210_T_C:211_L_N:211_L_CA -- 211_L_N:211_L_CA:211_L_C ] 
        [ [ 0.9440186246265643, 0.30880214609802714, -0.11606063468473155, -11.448370930248334 ], [ -0.3097717564868382, 0.9507578909206633, 0.010044487799343783, 3.5573829489349564 ], [ 0.11344732364066175, 0.026470123107851722, 0.9931913397434686, 20.994651778885675 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 210, ' ') T sidechain

     ],
     [ // 58 : (' ', 211, ' ') L backbone
      [ 141.08659, 232, 233, 0, 0, 1,     // 211_L_N:211_L_CA:211_L_C:211_L_O [ 211_L_N:211_L_CA:211_L_C -- 211_L_CA:211_L_C:211_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -37.96910, 232, 234, 0, 0, 1,     // 211_L_N:211_L_CA:211_L_C:212_E_N [ 211_L_N:211_L_CA:211_L_C -- 211_L_CA:211_L_C:212_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.83876, 234, 235, 0, 0, 1,     // 211_L_CA:211_L_C:212_E_N:212_E_CA [ 211_L_CA:211_L_C:212_E_N -- 211_L_C:212_E_N:212_E_CA ] 
        [ [ 0.34820499853902975, 0.6152364301385927, 0.7072746383285312, 0.0 ], [ -0.2717452874674142, 0.7883426507745984, -0.551969531504169, 0.0 ], [ -0.8971665273134563, 0.0, 0.4416924521296615, 15.265741942870875 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.87931, 235, 236, 0, 0, 1,     // 211_L_C:212_E_N:212_E_CA:212_E_C [ 211_L_C:212_E_N:212_E_CA -- 212_E_N:212_E_CA:212_E_C ] 
        [ [ -0.7827147633975121, -0.6162139219677485, 0.08739566083451894, 9.423926472466876 ], [ 0.612760800308706, -0.7875747758799447, -0.06519335857807085, -7.354597490205867 ], [ 0.10900367317026871, 0.0025248308419441806, 0.9940381403470436, 21.150976569904174 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 211, ' ') L sidechain

     ],
     [ // 59 : (' ', 212, ' ') E backbone
      [ 133.73077, 236, 237, 0, 0, 1,     // 212_E_N:212_E_CA:212_E_C:212_E_O [ 212_E_N:212_E_CA:212_E_C -- 212_E_CA:212_E_C:212_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -44.96631, 236, 238, 0, 0, 1,     // 212_E_N:212_E_CA:212_E_C:213_E_N [ 212_E_N:212_E_CA:212_E_C -- 212_E_CA:212_E_C:213_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.03835, 238, 239, 0, 0, 1,     // 212_E_CA:212_E_C:213_E_N:213_E_CA [ 212_E_CA:212_E_C:213_E_N -- 212_E_C:213_E_N:213_E_CA ] 
        [ [ 0.3142791648801956, 0.7066908898448744, 0.6338900478256986, 0.0 ], [ -0.3139097983047052, 0.7075224280616549, -0.6331450484036337, 0.0 ], [ -0.8959292634189968, 0.0, 0.44419675252020235, 15.289646887271747 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -65.77991, 239, 240, 0, 0, 1,     // 212_E_C:213_E_N:213_E_CA:213_E_C [ 212_E_C:213_E_N:213_E_CA -- 213_E_N:213_E_CA:213_E_C ] 
        [ [ -0.6982258943658629, -0.7118659602178564, 0.07567995190339522, 8.448049754325726 ], [ 0.7098657576335566, -0.7021543682750916, -0.05540622033259162, -8.438120915395162 ], [ 0.09258081105892023, 0.015036548650418929, 0.9955916309553612, 21.209595116262456 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 212, ' ') E sidechain

     ],
     [ // 60 : (' ', 213, ' ') E backbone
      [ 144.29057, 240, 241, 0, 0, 1,     // 213_E_N:213_E_CA:213_E_C:213_E_O [ 213_E_N:213_E_CA:213_E_C -- 213_E_CA:213_E_C:213_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  74.43916, 240, 242, 0, 0, 1,     // 213_E_N:213_E_CA:213_E_C:216_T_N [ 213_E_N:213_E_CA:213_E_C -- 213_E_CA:213_E_C:216_T_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  73.89487, 242, 243, 0, 0, 1,     // 213_E_CA:213_E_C:216_T_N:216_T_CA [ 213_E_CA:213_E_C:216_T_N -- 213_E_C:216_T_N:216_T_CA ] 
        [ [ 0.07534939345806434, -0.9633461399869134, 0.25746200394974017, 0.0 ], [ 0.2705850683764946, 0.2682614668048958, 0.924564495424399, 0.0 ], [ -0.9597427726621283, 0.0, 0.28088041997051066, 15.289168263607339 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  77.41879, 243, 244, 0, 0, 1,     // 213_E_C:216_T_N:216_T_CA:216_T_C [ 213_E_C:216_T_N:216_T_CA -- 216_T_N:216_T_CA:216_T_C ] 
        [ [ -0.5354137977144197, -0.33962519672575897, -0.7732960564790254, 9.782291239433013 ], [ -0.7678680721573918, -0.1855498352572264, 0.6131475209093703, 35.12890844136548 ], [ -0.3517253032955023, 0.9220769948247304, -0.16144140310552693, 25.961244257167024 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 213, ' ') E sidechain

     ],
     [ // 61 : (' ', 216, ' ') T backbone
      [ 135.89665, 244, 245, 0, 0, 1,     // 216_T_N:216_T_CA:216_T_C:216_T_O [ 216_T_N:216_T_CA:216_T_C -- 216_T_CA:216_T_C:216_T_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -41.95362, 244, 246, 0, 0, 1,     // 216_T_N:216_T_CA:216_T_C:217_A_N [ 216_T_N:216_T_CA:216_T_C -- 216_T_CA:216_T_C:217_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.47121, 246, 247, 0, 0, 1,     // 216_T_CA:216_T_C:217_A_N:217_A_CA [ 216_T_CA:216_T_C:217_A_N -- 216_T_C:217_A_N:217_A_CA ] 
        [ [ 0.32410803605027355, 0.668528853299139, 0.6693453169136039, 0.0 ], [ -0.29135349300874447, 0.7436862055373478, -0.6017008972945488, 0.0 ], [ -0.9000372898270593, 0.0, 0.4358128921002245, 15.274324612154437 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.28170, 247, 248, 0, 0, 1,     // 216_T_C:217_A_N:217_A_CA:217_A_C [ 216_T_C:217_A_N:217_A_CA -- 217_A_N:217_A_CA:217_A_C ] 
        [ [ -0.742117636146907, -0.665509180550481, 0.07974299168425013, 8.900958662310343 ], [ 0.6604701133470079, -0.7463434418639011, -0.0821626202243827, -8.001422701497907 ], [ 0.114195636935557, -0.008306466744228734, 0.9934235547414353, 21.06976747111556 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 216, ' ') T sidechain

     ],
     [ // 62 : (' ', 217, ' ') A backbone
      [ 138.19832, 248, 249, 0, 0, 1,     // 217_A_N:217_A_CA:217_A_C:217_A_O [ 217_A_N:217_A_CA:217_A_C -- 217_A_CA:217_A_C:217_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -39.92446, 248, 250, 0, 0, 1,     // 217_A_N:217_A_CA:217_A_C:218_C_N [ 217_A_N:217_A_CA:217_A_C -- 217_A_CA:217_A_C:218_C_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -178.73446, 250, 251, 0, 0, 1,     // 217_A_CA:217_A_C:218_C_N:218_C_CA [ 217_A_CA:217_A_C:218_C_N -- 217_A_C:218_C_N:218_C_CA ] 
        [ [ 0.3498419522490053, 0.6417771329554465, 0.6824461297876132, 0.0 ], [ -0.2927671686635241, 0.7668911993337041, -0.5711088104181693, 0.0 ], [ -0.8898865059092358, 0.0, 0.45618198846584423, 15.284906249934286 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -105.98237, 251, 252, 0, 0, 1,     // 217_A_C:218_C_N:218_C_CA:218_C_C [ 217_A_C:218_C_N:218_C_CA -- 218_C_N:218_C_CA:218_C_C ] 
        [ [ -0.7718850530028426, -0.6338939802697146, 0.048702019759119676, 9.082533100211847 ], [ 0.6308315040599266, -0.7731701888880549, -0.06526463437651538, -7.6007679552083305 ], [ 0.07902580867216531, -0.019654027389525975, 0.9966788052181516, 21.356137119469643 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 217, ' ') A sidechain

     ],
     [ // 63 : (' ', 218, ' ') C backbone
      [ -161.79296, 252, 253, 0, 0, 1,     // 218_C_N:218_C_CA:218_C_C:218_C_O [ 218_C_N:218_C_CA:218_C_C -- 218_C_CA:218_C_C:218_C_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  19.91767, 252, 254, 0, 0, 1,     // 218_C_N:218_C_CA:218_C_C:219_Q_N [ 218_C_N:218_C_CA:218_C_C -- 218_C_CA:218_C_C:219_Q_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.27934, 254, 255, 0, 0, 1,     // 218_C_CA:218_C_C:219_Q_N:219_Q_CA [ 218_C_CA:218_C_C:219_Q_N -- 218_C_C:219_Q_N:219_Q_CA ] 
        [ [ 0.43442209252845654, -0.3406694974449509, 0.8337995796555678, 0.0 ], [ 0.15741013985798052, 0.9401831170099815, 0.3021220851895578, 0.0 ], [ -0.8868480666907315, 0.0, 0.4620611502895606, 15.243517608489357 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.05711, 255, 256, 0, 0, 1,     // 218_C_C:219_Q_N:219_Q_CA:219_Q_C [ 218_C_C:219_Q_N:219_Q_CA -- 219_Q_N:219_Q_CA:219_Q_C ] 
        [ [ -0.936218610214449, 0.3461065290617864, 0.06086858326697517, 11.09185208217585 ], [ -0.34612428633759607, -0.9381289089639229, 0.010589078026421842, 4.0190635273112045 ], [ 0.0607675266521164, -0.011154403030313306, 0.9980896187204943, 21.390215276589586 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 218, ' ') C sidechain

     ],
     [ // 64 : (' ', 219, ' ') Q backbone
      [ -30.68307, 256, 257, 0, 0, 1,     // 219_Q_N:219_Q_CA:219_Q_C:219_Q_O [ 219_Q_N:219_Q_CA:219_Q_C -- 219_Q_CA:219_Q_C:219_Q_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 148.55963, 256, 258, 0, 0, 1,     // 219_Q_N:219_Q_CA:219_Q_C:220_G_N [ 219_Q_N:219_Q_CA:219_Q_C -- 219_Q_CA:219_Q_C:220_G_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.98219, 258, 259, 0, 0, 1,     // 219_Q_CA:219_Q_C:220_G_N:220_G_CA [ 219_Q_CA:219_Q_C:220_G_N -- 219_Q_C:220_G_N:220_G_CA ] 
        [ [ -0.37368104208034436, -0.5216109279120162, -0.7669970786596996, 0.0 ], [ 0.2284574433276205, -0.8531834737515521, 0.46891913663812196, 0.0 ], [ -0.8989825779056873, 0.0, 0.4379843885597351, 15.31298354339182 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 152.93178, 259, 260, 0, 0, 1,     // 219_Q_C:220_G_N:220_G_CA:220_G_C [ 219_Q_C:220_G_N:220_G_CA -- 220_G_N:220_G_CA:220_G_C ] 
        [ [ 0.849868695193912, 0.5217270431291521, -0.0743242450154019, -10.208783173799542 ], [ -0.5197717926560949, 0.8531124277612984, 0.04512725516019098, 6.24134553465203 ], [ 0.08695104650587629, 0.0002794046087920323, 0.9962125463196082, 21.14258548154903 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 219, ' ') Q sidechain

     ],
     [ // 65 : (' ', 220, ' ') G backbone
      [   5.90138, 260, 261, 0, 0, 1,     // 220_G_N:220_G_CA:220_G_C:220_G_O [ 220_G_N:220_G_CA:220_G_C -- 220_G_CA:220_G_C:220_G_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 220, ' ') G sidechain
      [ -174.49637, 260, 262, 0, 0, 1,     // 220_G_N:220_G_CA:220_G_C:220_G_OXT [ 220_G_N:220_G_CA:220_G_C -- 220_G_CA:220_G_C:220_G_OXT ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ]   ],
  ],
   [  //hedra
     [  14.54228, 113.45223,  15.34123, "Nsb", "Csb", "Cdb", 1, 1, 1, StdBond, StdBond, "D", 152, "NCAC" ], // (152_D_N, 152_D_CA, 152_D_C)
     [  15.34123, 120.62015,  12.32345, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 152, "CACO" ], // (152_D_CA, 152_D_C, 152_D_O)
     [  15.34123, 116.22613,  13.31165, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 152, "CACN" ], // (152_D_CA, 152_D_C, 153_I_N)
     [  13.31165, 122.90262,  14.65681, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 152, "CNCA" ], // (152_D_C, 153_I_N, 153_I_CA)
     [  14.65681, 114.34625,  15.34136, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "I", 153, "NCAC" ], // (153_I_N, 153_I_CA, 153_I_C)
     [  15.34136, 120.69603,  12.29753, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "I", 153, "CACO" ], // (153_I_CA, 153_I_C, 153_I_O)
     [  15.34136, 115.94010,  13.32844, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "I", 153, "CACN" ], // (153_I_CA, 153_I_C, 154_R_N)
     [  13.32844, 123.76076,  14.64516, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "I", 153, "CNCA" ], // (153_I_C, 154_R_N, 154_R_CA)
     [  14.64516, 107.28950,  15.23561, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "R", 154, "NCAC" ], // (154_R_N, 154_R_CA, 154_R_C)
     [  15.23561, 120.28015,  12.28659, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "R", 154, "CACO" ], // (154_R_CA, 154_R_C, 154_R_O)
     [  15.23561, 116.97304,  13.30194, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "R", 154, "CACN" ], // (154_R_CA, 154_R_C, 155_Q_N)
     [  13.30194, 121.41726,  14.62278, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "R", 154, "CNCA" ], // (154_R_C, 155_Q_N, 155_Q_CA)
     [  14.62278, 113.90657,  15.20362, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Q", 155, "NCAC" ], // (155_Q_N, 155_Q_CA, 155_Q_C)
     [  15.20362, 120.57797,  12.31903, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Q", 155, "CACO" ], // (155_Q_CA, 155_Q_C, 155_Q_O)
     [  15.20362, 116.39775,  13.29790, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Q", 155, "CACN" ], // (155_Q_CA, 155_Q_C, 156_G_N)
     [  13.29790, 120.28909,  14.51448, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Q", 155, "CNCA" ], // (155_Q_C, 156_G_N, 156_G_CA)
     [  14.51448, 109.52034,  15.22913, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "G", 156, "NCAC" ], // (156_G_N, 156_G_CA, 156_G_C)
     [  15.22913, 119.64706,  12.34741, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "G", 156, "CACO" ], // (156_G_CA, 156_G_C, 156_G_O)
     [  15.22913, 118.74633,  13.44059, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "G", 156, "CACN" ], // (156_G_CA, 156_G_C, 157_P_N)
     [  13.44059, 121.82570,  14.69110, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "G", 156, "CNCA" ], // (156_G_C, 157_P_N, 157_P_CA)
     [  14.69110, 114.43787,  15.32228, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "P", 157, "NCAC" ], // (157_P_N, 157_P_CA, 157_P_C)
     [  15.32228, 120.02921,  12.30497, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "P", 157, "CACO" ], // (157_P_CA, 157_P_C, 157_P_O)
     [  15.32228, 117.09641,  13.30775, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "P", 157, "CACN" ], // (157_P_CA, 157_P_C, 158_K_N)
     [  13.30775, 121.97739,  14.63850, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "P", 157, "CNCA" ], // (157_P_C, 158_K_N, 158_K_CA)
     [  14.63850, 112.55477,  15.29919, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 158, "NCAC" ], // (158_K_N, 158_K_CA, 158_K_C)
     [  15.29919, 120.71090,  12.33452, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 158, "CACO" ], // (158_K_CA, 158_K_C, 158_K_O)
     [  15.29919, 116.33389,  13.30991, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 158, "CACN" ], // (158_K_CA, 158_K_C, 159_E_N)
     [  13.30991, 122.06192,  14.59160, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 158, "CNCA" ], // (158_K_C, 159_E_N, 159_E_CA)
     [  14.59160, 110.80737,  15.26349, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 159, "NCAC" ], // (159_E_N, 159_E_CA, 159_E_C)
     [  15.26349, 120.32213,  12.28516, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 159, "CACO" ], // (159_E_CA, 159_E_C, 159_E_O)
     [  15.26349, 117.64322,  13.47035, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 159, "CACN" ], // (159_E_CA, 159_E_C, 160_P_N)
     [  13.47035, 121.20338,  14.68221, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 159, "CNCA" ], // (159_E_C, 160_P_N, 160_P_CA)
     [  14.68221, 109.00634,  15.28848, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "P", 160, "NCAC" ], // (160_P_N, 160_P_CA, 160_P_C)
     [  15.28848, 120.10565,  12.33209, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "P", 160, "CACO" ], // (160_P_CA, 160_P_C, 160_P_O)
     [  15.28848, 116.83124,  13.30082, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "P", 160, "CACN" ], // (160_P_CA, 160_P_C, 161_F_N)
     [  13.30082, 121.67711,  14.61247, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "P", 160, "CNCA" ], // (160_P_C, 161_F_N, 161_F_CA)
     [  14.61247, 114.21495,  15.31893, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "F", 161, "NCAC" ], // (161_F_N, 161_F_CA, 161_F_C)
     [  15.31893, 120.12146,  12.29607, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "F", 161, "CACO" ], // (161_F_CA, 161_F_C, 161_F_O)
     [  15.31893, 116.91942,  13.31399, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "F", 161, "CACN" ], // (161_F_CA, 161_F_C, 162_R_N)
     [  13.31399, 121.87774,  14.61196, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "F", 161, "CNCA" ], // (161_F_C, 162_R_N, 162_R_CA)
     [  14.61196, 112.27930,  15.35244, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "R", 162, "NCAC" ], // (162_R_N, 162_R_CA, 162_R_C)
     [  15.35244, 120.41181,  12.29248, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "R", 162, "CACO" ], // (162_R_CA, 162_R_C, 162_R_O)
     [  15.35244, 117.01295,  13.31699, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "R", 162, "CACN" ], // (162_R_CA, 162_R_C, 163_D_N)
     [  13.31699, 121.99159,  14.66164, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "R", 162, "CNCA" ], // (162_R_C, 163_D_N, 163_D_CA)
     [  14.66164, 110.65292,  15.29642, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "D", 163, "NCAC" ], // (163_D_N, 163_D_CA, 163_D_C)
     [  15.29642, 121.11136,  12.31612, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 163, "CACO" ], // (163_D_CA, 163_D_C, 163_D_O)
     [  15.29642, 115.53934,  13.27650, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 163, "CACN" ], // (163_D_CA, 163_D_C, 164_Y_N)
     [  13.27650, 122.85088,  14.54053, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 163, "CNCA" ], // (163_D_C, 164_Y_N, 164_Y_CA)
     [  14.54053, 108.55839,  15.31544, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Y", 164, "NCAC" ], // (164_Y_N, 164_Y_CA, 164_Y_C)
     [  15.31544, 120.55576,  12.29710, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Y", 164, "CACO" ], // (164_Y_CA, 164_Y_C, 164_Y_O)
     [  15.31544, 116.40067,  13.32446, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Y", 164, "CACN" ], // (164_Y_CA, 164_Y_C, 165_V_N)
     [  13.32446, 122.33029,  14.63621, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Y", 164, "CNCA" ], // (164_Y_C, 165_V_N, 165_V_CA)
     [  14.63621, 110.52887,  15.25646, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 165, "NCAC" ], // (165_V_N, 165_V_CA, 165_V_C)
     [  15.25646, 120.40551,  12.30967, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 165, "CACO" ], // (165_V_CA, 165_V_C, 165_V_O)
     [  15.25646, 116.59498,  13.31125, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 165, "CACN" ], // (165_V_CA, 165_V_C, 166_D_N)
     [  13.31125, 121.09761,  14.58892, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 165, "CNCA" ], // (165_V_C, 166_D_N, 166_D_CA)
     [  14.58892, 111.97087,  15.30166, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "D", 166, "NCAC" ], // (166_D_N, 166_D_CA, 166_D_C)
     [  15.30166, 120.20792,  12.24777, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 166, "CACO" ], // (166_D_CA, 166_D_C, 166_D_O)
     [  15.30166, 117.54498,  13.29855, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 166, "CACN" ], // (166_D_CA, 166_D_C, 167_R_N)
     [  13.29855, 120.60145,  14.63877, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 166, "CNCA" ], // (166_D_C, 167_R_N, 167_R_CA)
     [  14.63877, 109.38437,  15.25910, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "R", 167, "NCAC" ], // (167_R_N, 167_R_CA, 167_R_C)
     [  15.25910, 121.09128,  12.28527, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "R", 167, "CACO" ], // (167_R_CA, 167_R_C, 167_R_O)
     [  15.25910, 115.49619,  13.28353, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "R", 167, "CACN" ], // (167_R_CA, 167_R_C, 168_F_N)
     [  13.28353, 123.54740,  14.58920, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "R", 167, "CNCA" ], // (167_R_C, 168_F_N, 168_F_CA)
     [  14.58920, 111.01423,  15.25376, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "F", 168, "NCAC" ], // (168_F_N, 168_F_CA, 168_F_C)
     [  15.25376, 120.56250,  12.32539, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "F", 168, "CACO" ], // (168_F_CA, 168_F_C, 168_F_O)
     [  15.25376, 116.33157,  13.30793, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "F", 168, "CACN" ], // (168_F_CA, 168_F_C, 169_Y_N)
     [  13.30793, 121.77860,  14.61639, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "F", 168, "CNCA" ], // (168_F_C, 169_Y_N, 169_Y_CA)
     [  14.61639, 114.49506,  15.19299, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Y", 169, "NCAC" ], // (169_Y_N, 169_Y_CA, 169_Y_C)
     [  15.19299, 120.69058,  12.29097, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Y", 169, "CACO" ], // (169_Y_CA, 169_Y_C, 169_Y_O)
     [  15.19299, 116.26569,  13.26264, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Y", 169, "CACN" ], // (169_Y_CA, 169_Y_C, 170_K_N)
     [  13.26264, 121.49445,  14.54154, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Y", 169, "CNCA" ], // (169_Y_C, 170_K_N, 170_K_CA)
     [  14.54154, 106.00104,  15.30124, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 170, "NCAC" ], // (170_K_N, 170_K_CA, 170_K_C)
     [  15.30124, 120.72046,  12.25891, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 170, "CACO" ], // (170_K_CA, 170_K_C, 170_K_O)
     [  15.30124, 116.24835,  13.32164, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 170, "CACN" ], // (170_K_CA, 170_K_C, 171_T_N)
     [  13.32164, 122.05530,  14.60394, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 170, "CNCA" ], // (170_K_C, 171_T_N, 171_T_CA)
     [  14.60394, 109.91574,  15.29160, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "T", 171, "NCAC" ], // (171_T_N, 171_T_CA, 171_T_C)
     [  15.29160, 121.00788,  12.32871, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "T", 171, "CACO" ], // (171_T_CA, 171_T_C, 171_T_O)
     [  15.29160, 116.15896,  13.28342, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "T", 171, "CACN" ], // (171_T_CA, 171_T_C, 172_L_N)
     [  13.28342, 121.60771,  14.52686, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "T", 171, "CNCA" ], // (171_T_C, 172_L_N, 172_L_CA)
     [  14.52686, 110.95722,  15.32350, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 172, "NCAC" ], // (172_L_N, 172_L_CA, 172_L_C)
     [  15.32350, 120.28555,  12.32092, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 172, "CACO" ], // (172_L_CA, 172_L_C, 172_L_O)
     [  15.32350, 116.91610,  13.33082, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 172, "CACN" ], // (172_L_CA, 172_L_C, 173_R_N)
     [  13.33082, 121.73838,  14.64622, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 172, "CNCA" ], // (172_L_C, 173_R_N, 173_R_CA)
     [  14.64622, 112.88582,  15.21073, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "R", 173, "NCAC" ], // (173_R_N, 173_R_CA, 173_R_C)
     [  15.21073, 120.48348,  12.33075, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "R", 173, "CACO" ], // (173_R_CA, 173_R_C, 173_R_O)
     [  15.21073, 116.28116,  13.31198, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "R", 173, "CACN" ], // (173_R_CA, 173_R_C, 174_A_N)
     [  13.31198, 122.45068,  14.59428, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "R", 173, "CNCA" ], // (173_R_C, 174_A_N, 174_A_CA)
     [  14.59428, 112.59010,  15.24342, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 174, "NCAC" ], // (174_A_N, 174_A_CA, 174_A_C)
     [  15.24342, 120.52952,  12.30084, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 174, "CACO" ], // (174_A_CA, 174_A_C, 174_A_O)
     [  15.24342, 116.49978,  13.29326, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 174, "CACN" ], // (174_A_CA, 174_A_C, 175_E_N)
     [  13.29326, 121.96845,  14.57514, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 174, "CNCA" ], // (174_A_C, 175_E_N, 175_E_CA)
     [  14.57514, 112.18276,  15.25529, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 175, "NCAC" ], // (175_E_N, 175_E_CA, 175_E_C)
     [  15.25529, 120.16285,  12.30751, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 175, "CACO" ], // (175_E_CA, 175_E_C, 175_E_O)
     [  15.25529, 116.67290,  13.33947, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 175, "CACN" ], // (175_E_CA, 175_E_C, 176_Q_N)
     [  13.33947, 122.79975,  14.65177, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 175, "CNCA" ], // (175_E_C, 176_Q_N, 176_Q_CA)
     [  14.65177, 110.85477,  15.22863, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Q", 176, "NCAC" ], // (176_Q_N, 176_Q_CA, 176_Q_C)
     [  15.22863, 121.23578,  12.33024, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Q", 176, "CACO" ], // (176_Q_CA, 176_Q_C, 176_Q_O)
     [  15.22863, 115.50816,  13.28330, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Q", 176, "CACN" ], // (176_Q_CA, 176_Q_C, 177_A_N)
     [  13.28330, 122.17376,  14.53619, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Q", 176, "CNCA" ], // (176_Q_C, 177_A_N, 177_A_CA)
     [  14.53619, 111.62180,  15.22839, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 177, "NCAC" ], // (177_A_N, 177_A_CA, 177_A_C)
     [  15.22839, 121.04479,  12.32033, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 177, "CACO" ], // (177_A_CA, 177_A_C, 177_A_O)
     [  15.22839, 115.87363,  13.27122, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 177, "CACN" ], // (177_A_CA, 177_A_C, 178_S_N)
     [  13.27122, 122.16103,  14.57784, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 177, "CNCA" ], // (177_A_C, 178_S_N, 178_S_CA)
     [  14.57784, 108.47536,  15.25646, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 178, "NCAC" ], // (178_S_N, 178_S_CA, 178_S_C)
     [  15.25646, 120.39027,  12.29977, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 178, "CACO" ], // (178_S_CA, 178_S_C, 178_S_O)
     [  15.25646, 116.94739,  13.28418, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 178, "CACN" ], // (178_S_CA, 178_S_C, 179_Q_N)
     [  13.28418, 120.97936,  14.59698, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 178, "CNCA" ], // (178_S_C, 179_Q_N, 179_Q_CA)
     [  14.59698, 110.92994,  15.28417, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Q", 179, "NCAC" ], // (179_Q_N, 179_Q_CA, 179_Q_C)
     [  15.28417, 120.54496,  12.32779, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Q", 179, "CACO" ], // (179_Q_CA, 179_Q_C, 179_Q_O)
     [  15.28417, 116.25597,  13.31214, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Q", 179, "CACN" ], // (179_Q_CA, 179_Q_C, 180_E_N)
     [  13.31214, 121.94458,  14.59266, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Q", 179, "CNCA" ], // (179_Q_C, 180_E_N, 180_E_CA)
     [  14.59266, 111.92690,  15.22120, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 180, "NCAC" ], // (180_E_N, 180_E_CA, 180_E_C)
     [  15.22120, 120.16493,  12.30421, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 180, "CACO" ], // (180_E_CA, 180_E_C, 180_E_O)
     [  15.22120, 116.88175,  13.30887, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 180, "CACN" ], // (180_E_CA, 180_E_C, 181_V_N)
     [  13.30887, 121.23235,  14.61270, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 180, "CNCA" ], // (180_E_C, 181_V_N, 181_V_CA)
     [  14.61270, 110.20384,  15.26821, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 181, "NCAC" ], // (181_V_N, 181_V_CA, 181_V_C)
     [  15.26821, 120.65420,  12.31038, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 181, "CACO" ], // (181_V_CA, 181_V_C, 181_V_O)
     [  15.26821, 116.53634,  13.28645, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 181, "CACN" ], // (181_V_CA, 181_V_C, 182_K_N)
     [  13.28645, 121.06826,  14.58320, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 181, "CNCA" ], // (181_V_C, 182_K_N, 182_K_CA)
     [  14.58320, 109.87763,  15.25626, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 182, "NCAC" ], // (182_K_N, 182_K_CA, 182_K_C)
     [  15.25626, 120.21760,  12.31697, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 182, "CACO" ], // (182_K_CA, 182_K_C, 182_K_O)
     [  15.25626, 116.56205,  13.33937, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 182, "CACN" ], // (182_K_CA, 182_K_C, 183_N_N)
     [  13.33937, 121.89817,  14.63992, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 182, "CNCA" ], // (182_K_C, 183_N_N, 183_N_CA)
     [  14.63992, 110.84403,  15.18493, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "N", 183, "NCAC" ], // (183_N_N, 183_N_CA, 183_N_C)
     [  15.18493, 120.10337,  12.26153, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "N", 183, "CACO" ], // (183_N_CA, 183_N_C, 183_N_O)
     [  15.18493, 116.65884,  13.31695, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "N", 183, "CACN" ], // (183_N_CA, 183_N_C, 184_W_N)
     [  13.31695, 121.85759,  14.59773, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "N", 183, "CNCA" ], // (183_N_C, 184_W_N, 184_W_CA)
     [  14.59773, 110.50903,  15.24341, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "W", 184, "NCAC" ], // (184_W_N, 184_W_CA, 184_W_C)
     [  15.24341, 120.50272,  12.30235, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "W", 184, "CACO" ], // (184_W_CA, 184_W_C, 184_W_O)
     [  15.24341, 127.42233,  32.68179, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "W", 184, "CACN" ], // (184_W_CA, 184_W_C, 186_T_N)
     [  32.68179, 139.18954,  14.66632, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "W", 184, "CNCA" ], // (184_W_C, 186_T_N, 186_T_CA)
     [  14.66632, 111.53921,  15.40964, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "T", 186, "NCAC" ], // (186_T_N, 186_T_CA, 186_T_C)
     [  15.40964, 120.25820,  12.37163, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "T", 186, "CACO" ], // (186_T_CA, 186_T_C, 186_T_O)
     [  15.40964, 116.39008,  13.37088, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "T", 186, "CACN" ], // (186_T_CA, 186_T_C, 187_E_N)
     [  13.37088, 123.38465,  14.70596, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "T", 186, "CNCA" ], // (186_T_C, 187_E_N, 187_E_CA)
     [  14.70596, 112.06115,  15.20403, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 187, "NCAC" ], // (187_E_N, 187_E_CA, 187_E_C)
     [  15.20403, 120.28115,  12.28975, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 187, "CACO" ], // (187_E_CA, 187_E_C, 187_E_O)
     [  15.20403, 116.42023,  13.31249, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 187, "CACN" ], // (187_E_CA, 187_E_C, 188_T_N)
     [  13.31249, 121.97345,  14.58094, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 187, "CNCA" ], // (187_E_C, 188_T_N, 188_T_CA)
     [  14.58094, 116.66723,  15.23286, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "T", 188, "NCAC" ], // (188_T_N, 188_T_CA, 188_T_C)
     [  15.23286, 120.60451,  12.31781, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "T", 188, "CACO" ], // (188_T_CA, 188_T_C, 188_T_O)
     [  15.23286, 116.13974,  13.30802, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "T", 188, "CACN" ], // (188_T_CA, 188_T_C, 189_L_N)
     [  13.30802, 122.53297,  14.54168, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "T", 188, "CNCA" ], // (188_T_C, 189_L_N, 189_L_CA)
     [  14.54168, 118.02114,  15.23827, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 189, "NCAC" ], // (189_L_N, 189_L_CA, 189_L_C)
     [  15.23827, 119.45463,  12.30366, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 189, "CACO" ], // (189_L_CA, 189_L_C, 189_L_O)
     [  15.23827, 118.32114,  13.31710, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 189, "CACN" ], // (189_L_CA, 189_L_C, 190_L_N)
     [  13.31710, 119.54300,  14.65072, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 189, "CNCA" ], // (189_L_C, 190_L_N, 190_L_CA)
     [  14.65072, 110.06257,  15.25307, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 190, "NCAC" ], // (190_L_N, 190_L_CA, 190_L_C)
     [  15.25307, 120.81610,  12.28155, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 190, "CACO" ], // (190_L_CA, 190_L_C, 190_L_O)
     [  15.25307, 116.03074,  13.29277, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 190, "CACN" ], // (190_L_CA, 190_L_C, 191_V_N)
     [  13.29277, 122.16139,  14.58032, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 190, "CNCA" ], // (190_L_C, 191_V_N, 191_V_CA)
     [  14.58032, 110.58199,  15.31034, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 191, "NCAC" ], // (191_V_N, 191_V_CA, 191_V_C)
     [  15.31034, 120.75243,  12.32039, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 191, "CACO" ], // (191_V_CA, 191_V_C, 191_V_O)
     [  15.31034, 116.05094,  13.30432, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 191, "CACN" ], // (191_V_CA, 191_V_C, 192_Q_N)
     [  13.30432, 122.75456,  14.61821, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 191, "CNCA" ], // (191_V_C, 192_Q_N, 192_Q_CA)
     [  14.61821, 114.85184,  15.25865, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Q", 192, "NCAC" ], // (192_Q_N, 192_Q_CA, 192_Q_C)
     [  15.25865, 120.23786,  12.28998, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Q", 192, "CACO" ], // (192_Q_CA, 192_Q_C, 192_Q_O)
     [  15.25865, 117.05245,  13.29317, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Q", 192, "CACN" ], // (192_Q_CA, 192_Q_C, 193_N_N)
     [  13.29317, 121.21041,  14.58403, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Q", 192, "CNCA" ], // (192_Q_C, 193_N_N, 193_N_CA)
     [  14.58403, 113.14353,  15.26659, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "N", 193, "NCAC" ], // (193_N_N, 193_N_CA, 193_N_C)
     [  15.26659, 120.47126,  12.28776, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "N", 193, "CACO" ], // (193_N_CA, 193_N_C, 193_N_O)
     [  15.26659, 116.53007,  13.30805, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "N", 193, "CACN" ], // (193_N_CA, 193_N_C, 194_A_N)
     [  13.30805, 122.15075,  14.61493, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "N", 193, "CNCA" ], // (193_N_C, 194_A_N, 194_A_CA)
     [  14.61493, 109.47068,  15.35186, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 194, "NCAC" ], // (194_A_N, 194_A_CA, 194_A_C)
     [  15.35186, 120.39739,  12.34513, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 194, "CACO" ], // (194_A_CA, 194_A_C, 194_A_O)
     [  15.35186, 117.26700,  13.32205, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 194, "CACN" ], // (194_A_CA, 194_A_C, 195_N_N)
     [  13.32205, 120.41949,  14.61776, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 194, "CNCA" ], // (194_A_C, 195_N_N, 195_N_CA)
     [  14.61776, 111.97136,  15.34245, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "N", 195, "NCAC" ], // (195_N_N, 195_N_CA, 195_N_C)
     [  15.34245, 120.36004,  12.28661, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "N", 195, "CACO" ], // (195_N_CA, 195_N_C, 195_N_O)
     [  15.34245, 117.76212,  13.45786, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "N", 195, "CACN" ], // (195_N_CA, 195_N_C, 196_P_N)
     [  13.45786, 123.39439,  14.69310, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "N", 195, "CNCA" ], // (195_N_C, 196_P_N, 196_P_CA)
     [  14.69310, 115.75464,  15.29755, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "P", 196, "NCAC" ], // (196_P_N, 196_P_CA, 196_P_C)
     [  15.29755, 120.35146,  12.31527, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "P", 196, "CACO" ], // (196_P_CA, 196_P_C, 196_P_O)
     [  15.29755, 116.70069,  13.28773, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "P", 196, "CACN" ], // (196_P_CA, 196_P_C, 197_D_N)
     [  13.28773, 121.80099,  14.59515, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "P", 196, "CNCA" ], // (196_P_C, 197_D_N, 197_D_CA)
     [  14.59515, 113.24398,  15.23121, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "D", 197, "NCAC" ], // (197_D_N, 197_D_CA, 197_D_C)
     [  15.23121, 120.75624,  12.31677, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 197, "CACO" ], // (197_D_CA, 197_D_C, 197_D_O)
     [  15.23121, 116.02402,  13.28568, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 197, "CACN" ], // (197_D_CA, 197_D_C, 198_C_N)
     [  13.28568, 121.88727,  14.58759, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 197, "CNCA" ], // (197_D_C, 198_C_N, 198_C_CA)
     [  14.58759, 115.35729,  15.31695, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "C", 198, "NCAC" ], // (198_C_N, 198_C_CA, 198_C_C)
     [  15.31695, 120.34473,  12.30296, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "C", 198, "CACO" ], // (198_C_CA, 198_C_C, 198_C_O)
     [  15.31695, 117.02594,  13.30291, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "C", 198, "CACN" ], // (198_C_CA, 198_C_C, 199_K_N)
     [  13.30291, 121.91423,  14.61692, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "C", 198, "CNCA" ], // (198_C_C, 199_K_N, 199_K_CA)
     [  14.61692, 113.82421,  15.27006, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 199, "NCAC" ], // (199_K_N, 199_K_CA, 199_K_C)
     [  15.27006, 120.36345,  12.33418, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 199, "CACO" ], // (199_K_CA, 199_K_C, 199_K_O)
     [  15.27006, 116.61228,  13.30630, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 199, "CACN" ], // (199_K_CA, 199_K_C, 200_T_N)
     [  13.30630, 121.72589,  14.63211, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 199, "CNCA" ], // (199_K_C, 200_T_N, 200_T_CA)
     [  14.63211, 110.16606,  15.28148, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "T", 200, "NCAC" ], // (200_T_N, 200_T_CA, 200_T_C)
     [  15.28148, 121.22192,  12.31772, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "T", 200, "CACO" ], // (200_T_CA, 200_T_C, 200_T_O)
     [  15.28148, 115.40689,  13.27671, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "T", 200, "CACN" ], // (200_T_CA, 200_T_C, 201_I_N)
     [  13.27671, 122.92136,  14.54422, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "T", 200, "CNCA" ], // (200_T_C, 201_I_N, 201_I_CA)
     [  14.54422, 111.15643,  15.27360, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "I", 201, "NCAC" ], // (201_I_N, 201_I_CA, 201_I_C)
     [  15.27360, 120.72744,  12.30728, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "I", 201, "CACO" ], // (201_I_CA, 201_I_C, 201_I_O)
     [  15.27360, 116.25105,  13.29794, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "I", 201, "CACN" ], // (201_I_CA, 201_I_C, 202_L_N)
     [  13.29794, 122.00866,  14.57238, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "I", 201, "CNCA" ], // (201_I_C, 202_L_N, 202_L_CA)
     [  14.57238, 114.47754,  15.23422, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 202, "NCAC" ], // (202_L_N, 202_L_CA, 202_L_C)
     [  15.23422, 120.31178,  12.31038, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 202, "CACO" ], // (202_L_CA, 202_L_C, 202_L_O)
     [  15.23422, 116.94946,  13.28438, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 202, "CACN" ], // (202_L_CA, 202_L_C, 203_K_N)
     [  13.28438, 121.07374,  14.60089, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 202, "CNCA" ], // (202_L_C, 203_K_N, 203_K_CA)
     [  14.60089, 110.00642,  15.28429, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 203, "NCAC" ], // (203_K_N, 203_K_CA, 203_K_C)
     [  15.28429, 120.81296,  12.34844, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 203, "CACO" ], // (203_K_CA, 203_K_C, 203_K_O)
     [  15.28429, 115.65182,  13.32521, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 203, "CACN" ], // (203_K_CA, 203_K_C, 204_A_N)
     [  13.32521, 123.07176,  14.60556, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 203, "CNCA" ], // (203_K_C, 204_A_N, 204_A_CA)
     [  14.60556, 114.22698,  15.35056, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 204, "NCAC" ], // (204_A_N, 204_A_CA, 204_A_C)
     [  15.35056, 120.43578,  12.34845, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 204, "CACO" ], // (204_A_CA, 204_A_C, 204_A_O)
     [  15.35056, 116.99568,  13.31111, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 204, "CACN" ], // (204_A_CA, 204_A_C, 205_L_N)
     [  13.31111, 121.42652,  14.60748, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 204, "CNCA" ], // (204_A_C, 205_L_N, 205_L_CA)
     [  14.60748, 110.18352,  15.26584, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 205, "NCAC" ], // (205_L_N, 205_L_CA, 205_L_C)
     [  15.26584, 120.72708,  12.31243, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 205, "CACO" ], // (205_L_CA, 205_L_C, 205_L_O)
     [  15.26584, 116.32498,  13.31176, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 205, "CACN" ], // (205_L_CA, 205_L_C, 206_G_N)
     [  13.31176, 120.73645,  14.52692, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 205, "CNCA" ], // (205_L_C, 206_G_N, 206_G_CA)
     [  14.52692, 112.03008,  15.14729, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "G", 206, "NCAC" ], // (206_G_N, 206_G_CA, 206_G_C)
     [  15.14729, 120.08096,  12.31230, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "G", 206, "CACO" ], // (206_G_CA, 206_G_C, 206_G_O)
     [  15.14729, 118.26019,  13.42183, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "G", 206, "CACN" ], // (206_G_CA, 206_G_C, 207_P_N)
     [  13.42183, 122.23208,  14.68762, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "G", 206, "CNCA" ], // (206_G_C, 207_P_N, 207_P_CA)
     [  14.68762, 115.83759,  15.31826, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "P", 207, "NCAC" ], // (207_P_N, 207_P_CA, 207_P_C)
     [  15.31826, 121.43279,  12.32966, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "P", 207, "CACO" ], // (207_P_CA, 207_P_C, 207_P_O)
     [  15.31826, 115.32328,  13.29650, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "P", 207, "CACN" ], // (207_P_CA, 207_P_C, 208_G_N)
     [  13.29650, 121.08552,  14.49437, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "P", 207, "CNCA" ], // (207_P_C, 208_G_N, 208_G_CA)
     [  14.49437, 116.55529,  15.16163, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "G", 208, "NCAC" ], // (208_G_N, 208_G_CA, 208_G_C)
     [  15.16163, 119.30449,  12.31274, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "G", 208, "CACO" ], // (208_G_CA, 208_G_C, 208_G_O)
     [  15.16163, 117.83217,  13.31859, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "G", 208, "CACN" ], // (208_G_CA, 208_G_C, 209_A_N)
     [  13.31859, 121.11717,  14.65673, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "G", 208, "CNCA" ], // (208_G_C, 209_A_N, 209_A_CA)
     [  14.65673, 108.75548,  15.26642, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 209, "NCAC" ], // (209_A_N, 209_A_CA, 209_A_C)
     [  15.26642, 120.56949,  12.33975, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 209, "CACO" ], // (209_A_CA, 209_A_C, 209_A_O)
     [  15.26642, 116.52060,  13.29321, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 209, "CACN" ], // (209_A_CA, 209_A_C, 210_T_N)
     [  13.29321, 121.08399,  14.56646, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 209, "CNCA" ], // (209_A_C, 210_T_N, 210_T_CA)
     [  14.56646, 112.10139,  15.21185, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "T", 210, "NCAC" ], // (210_T_N, 210_T_CA, 210_T_C)
     [  15.21185, 120.68763,  12.30813, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "T", 210, "CACO" ], // (210_T_CA, 210_T_C, 210_T_O)
     [  15.21185, 115.75124,  13.31018, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "T", 210, "CACN" ], // (210_T_CA, 210_T_C, 211_L_N)
     [  13.31018, 122.27730,  14.56928, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "T", 210, "CNCA" ], // (210_T_C, 211_L_N, 211_L_CA)
     [  14.56928, 108.43813,  15.26574, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 211, "NCAC" ], // (211_L_N, 211_L_CA, 211_L_C)
     [  15.26574, 120.53873,  12.29074, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 211, "CACO" ], // (211_L_CA, 211_L_C, 211_L_O)
     [  15.26574, 116.21192,  13.32428, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 211, "CACN" ], // (211_L_CA, 211_L_C, 212_E_N)
     [  13.32428, 122.46991,  14.59793, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 211, "CNCA" ], // (211_L_C, 212_E_N, 212_E_CA)
     [  14.59793, 111.14993,  15.28965, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 212, "NCAC" ], // (212_E_N, 212_E_CA, 212_E_C)
     [  15.28965, 120.55874,  12.32823, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 212, "CACO" ], // (212_E_CA, 212_E_C, 212_E_O)
     [  15.28965, 116.37196,  13.32731, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 212, "CACN" ], // (212_E_CA, 212_E_C, 213_E_N)
     [  13.32731, 121.68787,  14.62456, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 212, "CNCA" ], // (212_E_C, 213_E_N, 213_E_CA)
     [  14.62456, 109.98364,  15.28917, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 213, "NCAC" ], // (213_E_N, 213_E_CA, 213_E_C)
     [  15.28917, 120.45960,  12.32717, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 213, "CACO" ], // (213_E_CA, 213_E_C, 213_E_O)
     [  15.28917, 106.31276,  37.99509, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 213, "CACN" ], // (213_E_CA, 213_E_C, 216_T_N)
     [  37.99509, 108.81142,  14.58733, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 213, "CNCA" ], // (213_E_C, 216_T_N, 216_T_CA)
     [  14.58733, 109.49424,  15.27432, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "T", 216, "NCAC" ], // (216_T_N, 216_T_CA, 216_T_C)
     [  15.27432, 120.84739,  12.31874, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "T", 216, "CACO" ], // (216_T_CA, 216_T_C, 216_T_O)
     [  15.27432, 115.83703,  13.29801, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "T", 216, "CACN" ], // (216_T_CA, 216_T_C, 217_A_N)
     [  13.29801, 122.39545,  14.57324, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "T", 216, "CNCA" ], // (216_T_C, 217_A_N, 217_A_CA)
     [  14.57324, 114.07797,  15.28491, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 217, "NCAC" ], // (217_A_N, 217_A_CA, 217_A_C)
     [  15.28491, 120.19127,  12.28033, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 217, "CACO" ], // (217_A_CA, 217_A_C, 217_A_O)
     [  15.28491, 117.14101,  13.30879, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 217, "CACN" ], // (217_A_CA, 217_A_C, 218_C_N)
     [  13.30879, 121.68013,  14.63380, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 217, "CNCA" ], // (217_A_C, 218_C_N, 218_C_CA)
     [  14.63380, 115.95453,  15.24352, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "C", 218, "NCAC" ], // (218_C_N, 218_C_CA, 218_C_C)
     [  15.24352, 119.76113,  12.31508, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "C", 218, "CACO" ], // (218_C_CA, 218_C_C, 218_C_O)
     [  15.24352, 117.52019,  13.30278, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "C", 218, "CACN" ], // (218_C_CA, 218_C_C, 219_Q_N)
     [  13.30278, 121.00613,  14.62498, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "C", 218, "CNCA" ], // (218_C_C, 219_Q_N, 219_Q_CA)
     [  14.62498, 112.89591,  15.31298, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Q", 219, "NCAC" ], // (219_Q_N, 219_Q_CA, 219_Q_C)
     [  15.31298, 120.96332,  12.31152, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Q", 219, "CACO" ], // (219_Q_CA, 219_Q_C, 219_Q_O)
     [  15.31298, 115.97535,  13.31007, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Q", 219, "CACN" ], // (219_Q_CA, 219_Q_C, 220_G_N)
     [  13.31007, 120.96358,  14.55115, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Q", 219, "CNCA" ], // (219_Q_C, 220_G_N, 220_G_CA)
     [  14.55115, 114.14807,  15.22995, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "G", 220, "NCAC" ], // (220_G_N, 220_G_CA, 220_G_C)
     [  15.22995, 119.13228,  12.51227, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "G", 220, "CACO" ], // (220_G_CA, 220_G_C, 220_G_O)
     [  15.22995, 118.44799,  12.49952, "Csb", "Cdb", "Osb", 0, 0, 1, 0, StdBond, "G", 220, "CACOXT" ], // (220_G_CA, 220_G_C, 220_G_OXT)
   ],

[  // chain - world transform for each residue
     [ 0, "152D", //(152_D_N, 152_D_CA, 152_D_C)
[ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 1, "153I", //(153_I_N, 153_I_CA, 153_I_C)
[ [ -0.09925244296660844, -0.9949352024582806, 0.015902687648962906, 12.105682919280593 ], [ 0.4714176900638687, -0.03294124480521226, 0.8812946362514238, -6.303102364028917 ], [ -0.8763072030173017, 0.09496745389860368, 0.47229955392738854, 35.779275354548616 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 2, "154R", //(154_R_N, 154_R_CA, 154_R_C)
[ [ 0.02032516563866696, 0.9804184507995578, -0.19587380879932825, 1.8559346259091143 ], [ -0.9558618966071663, -0.03838129239002528, -0.2912986629029091, 21.104234037407192 ], [ -0.29311247373032767, 0.19314900394841916, 0.9363645337246715, 60.52103223632049 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 3, "155Q", //(155_Q_N, 155_Q_CA, 155_Q_C)
[ [ -0.2567371036356283, -0.34722394527944744, -0.9019543178240037, 3.6171425975112186 ], [ 0.37952110595437466, -0.8944917852861113, 0.23632218725890639, 19.128624806842826 ], [ -0.8888474502140156, -0.2816380263392012, 0.3614280459065157, 98.55618453915693 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 4, "156G", //(156_G_N, 156_G_CA, 156_G_C)
[ [ 0.988678655642436, -0.10870949831582473, 0.10342514613471546, -29.91702162867152 ], [ 0.024859236292028356, 0.7984250164721195, 0.6015808436465317, 15.46536391802208 ], [ -0.1479747757154666, -0.592199069610395, 0.7920882070227658, 115.74170820683358 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 5, "157P", //(157_P_N, 157_P_CA, 157_P_C)
[ [ 0.025589227606157044, 0.9819757026127993, 0.18726695092466084, -37.68669331972592 ], [ -0.8862986787978724, 0.10893571497832844, -0.4501196085079012, 41.94258375406096 ], [ -0.46240657801513396, -0.15445623807492992, 0.8731113486423564, 142.5549088589193 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 6, "158K", //(158_K_N, 158_K_CA, 158_K_C)
[ [ 0.06632629789121929, 0.7996839220853463, 0.596746551700337, -35.584307995432184 ], [ 0.24127637272543273, 0.5674614447458511, -0.7872567692255236, 14.490065038857995 ], [ -0.9681872412774831, 0.1961966704231982, -0.15530721922160495, 169.14622676283892 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 7, "159E", //(159_E_N, 159_E_CA, 159_E_C)
[ [ -0.30757177164963356, -0.227715032674088, 0.9238752454625813, -11.471079560464869 ], [ 0.8668719297661792, -0.4674070408327971, 0.17338891418827318, -9.091099707388027 ], [ 0.3923425323206319, 0.8542110524186438, 0.34116098144141127, 151.32164302027755 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 8, "160P", //(160_P_N, 160_P_CA, 160_P_C)
[ [ -0.7607123886964525, -0.1571359337797711, 0.6297816764553384, 22.072771715204308 ], [ -0.06408795980734479, 0.9836965957039848, 0.1680289885945968, -14.40915057918922 ], [ -0.6459174831907474, 0.08746031051597278, -0.7583807084789145, 169.19803546618598 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 9, "161F", //(161_F_N, 161_F_CA, 161_F_C)
[ [ -0.5058762647196147, -0.71225356055778, 0.48660463444781865, 50.97431717668728 ], [ -0.6096791771882328, -0.10383951157700932, -0.7858171904066031, 0.5331737992123156 ], [ 0.6102298793867846, -0.6941989781814815, -0.38171622050336973, 149.23649466379553 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 10, "162R", //(162_R_N, 162_R_CA, 162_R_C)
[ [ -0.7540523388236562, 0.3747264764540624, -0.5394303830513227, 69.82537656234571 ], [ 0.6534002148704969, 0.5115963625500581, -0.5579760936032454, -32.68651658671357 ], [ 0.06688220631649507, -0.7732071065826234, -0.6306167939472952, 147.59567167359796 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 11, "163D", //(163_D_N, 163_D_CA, 163_D_C)
[ [ 0.6268528662339278, 0.7760197825817529, -0.06963318990308079, 39.5322428175497 ], [ 0.7234344473910276, -0.5465225897124798, 0.42184791010742617, -49.99047026420948 ], [ 0.28930621220835906, -0.31481161982324213, -0.9039886943993926, 131.72271850355472 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 12, "164Y", //(164_Y_N, 164_Y_CA, 164_Y_C)
[ [ 0.33813315992816495, -0.37714113986681713, 0.8622241743171853, 37.015646527653956 ], [ -0.6015013506591789, -0.7912341612784004, -0.110202664129548, -22.503120610583103 ], [ 0.7237831797659972, -0.48136583036815594, -0.49439341221520955, 105.5429829686337 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 13, "165V", //(165_V_N, 165_V_CA, 165_V_C)
[ [ -0.9454747749563309, -0.31869087278564273, 0.06718316399516043, 74.64868111525043 ], [ -0.2080934676877374, 0.4324092525091614, -0.8773365073050141, -25.42983478155411 ], [ 0.2505485155153946, -0.8434799143711684, -0.4751495295435306, 99.42207460618378 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 14, "166D", //(166_D_N, 166_D_CA, 166_D_C)
[ [ -0.2486241157254675, 0.7977652674411239, -0.5493237908049508, 70.97264111634362 ], [ 0.9620585522906042, 0.2691913347405034, -0.044490080525133945, -62.657844146032 ], [ 0.11238056346288715, -0.539542957849714, -0.8344243558229504, 92.38284324894562 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 15, "167R", //(167_R_N, 167_R_CA, 167_R_C)
[ [ 0.774623524970798, 0.4381535102585997, 0.4560481290498654, 42.216774294411245 ], [ 0.17911135274335474, -0.8435816490941337, 0.5062500613630259, -57.35604628107819 ], [ 0.606529074225021, -0.3104698097395125, -0.7319364585536438, 67.80705059444692 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 16, "168F", //(168_F_N, 168_F_CA, 168_F_C)
[ [ -0.23634399997369984, -0.6121466496297729, 0.7545979015498755, 61.93910086145467 ], [ -0.8282943804188538, -0.2790946432486065, -0.48583392170213885, -29.402927792092918 ], [ 0.5080058395756675, -0.7398531337076143, -0.44107528552396164, 51.03433106640939 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 17, "169Y", //(169_Y_N, 169_Y_CA, 169_Y_C)
[ [ -0.969849514424169, -0.025626706869188796, -0.24235344286030777, 93.04274679092644 ], [ 0.19034661441394557, 0.5413450996182966, -0.8189710919809108, -51.080024067030934 ], [ 0.1521843807765847, -0.8404098732251408, -0.5201453251095796, 47.16918272034701 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 18, "170K", //(170_K_N, 170_K_CA, 170_K_C)
[ [ 0.30119620268688396, 0.9413430625366571, -0.15216466771596043, 72.6745587177187 ], [ 0.9226141625564811, -0.24736199758937394, 0.2959647769563086, -80.54856990269386 ], [ 0.2409646333742869, -0.2295327444237568, -0.9430009356835728, 34.73327604158992 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 19, "171T", //(171_T_N, 171_T_CA, 171_T_C)
[ [ 0.49935215214492634, 0.07707807117979701, 0.8629637298817561, 61.87508624349655 ], [ -0.31772811179284904, -0.9103501033049941, 0.2651632259366622, -58.50322777861079 ], [ 0.8060373906493304, -0.4065976639822082, -0.4300954132741131, 5.51726743300825 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 20, "172L", //(172_L_N, 172_L_CA, 172_L_C)
[ [ -0.7566021138167341, -0.38853564101808097, 0.5259213791306823, 97.49408187413518 ], [ -0.607365581767269, 0.11967587217258047, -0.7853564386342289, -45.67502509010753 ], [ 0.2421988675707565, -0.913628885969791, -0.32653019350290113, 2.100619809123718 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 21, "173R", //(173_R_N, 173_R_CA, 173_R_C)
[ [ -0.6475767339992765, 0.7074949359500392, -0.28301111142123064, 111.59347350382255 ], [ 0.7569826565577091, 0.5547451955883361, -0.345304250831773, -81.20037549785094 ], [ -0.08730195446645855, -0.43784550194869953, -0.8948014780774716, 0.15099680000213933 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 22, "174A", //(174_A_N, 174_A_CA, 174_A_C)
[ [ 0.5977028689114834, 0.6701617872154645, 0.44003915672489485, 88.90494592532463 ], [ 0.5175619852861897, -0.7417152764698137, 0.4266005626319845, -89.27675594024244 ], [ 0.6122751602683225, -0.027232840606192698, -0.7901756137156586, -29.35786475732561 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 23, "175E", //(175_E_N, 175_E_CA, 175_E_C)
[ [ -0.1059540969911183, -0.0952581867854883, 0.9897977607476882, 109.27107555821394 ], [ -0.5902783496966836, -0.7950191518467693, -0.1396997425774411, -64.37953730207113 ], [ 0.800215720421785, -0.5990579488223027, 0.028006691033732736, -49.733765334762474 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 24, "176Q", //(176_Q_N, 176_Q_CA, 176_Q_C)
[ [ -0.8366332562675416, -0.10751154384153427, 0.5371089856332512, 143.20846825130084 ], [ 0.5347483422154274, 0.05217182891499609, 0.843399259405364, -78.00784936174874 ], [ -0.1186971145606424, 0.9928340083863335, 0.013842932727695446, -38.47509661196203 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 25, "177A", //(177_A_N, 177_A_CA, 177_A_C)
[ [ -0.12599504975254044, 0.3642911696180483, 0.9227227054625724, 153.38989926171797 ], [ -0.9905612078107628, -0.09681021444106713, -0.09703749770247208, -42.68851076948795 ], [ 0.053979079450631025, -0.9262395619483055, 0.3730503087567058, -29.08690007674695 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 26, "178S", //(178_S_N, 178_S_CA, 178_S_C)
[ [ -0.984142029988627, -0.1697770958950502, 0.05138290104020278, 190.58175738684332 ], [ 0.1727880428160671, -0.8520506750927089, 0.4941193573761828, -36.5282596768807 ], [ -0.04010931400133012, 0.4951619783298207, 0.8678743331536202, -24.29557986841011 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 27, "179Q", //(179_Q_N, 179_Q_CA, 179_Q_C)
[ [ -0.21371818095940276, -0.9521318195691179, -0.21856243339468168, 205.0230072628012 ], [ -0.14798056229113854, -0.18959639627115998, 0.9706466708875006, -21.34520668371487 ], [ -0.9656222306427097, 0.23978783264584122, -0.10037680510759252, 7.497763748455297 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 28, "180E", //(180_E_N, 180_E_CA, 180_E_C)
[ [ -0.20977256497569816, -0.06839680557168072, -0.975354985618624, 204.38246304133057 ], [ -0.964571597204387, 0.17770735623979889, 0.19499161367822693, 13.744466703199754 ], [ 0.15999095240073477, 0.9817036072694502, -0.10325174392742674, -7.385703101064365 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 29, "181V", //(181_V_N, 181_V_CA, 181_V_C)
[ [ 0.870136405225635, 0.1070988822051808, -0.4810327075484683, 168.0170476851748 ], [ -0.1599952027151097, 0.9846188101311767, -0.07019497021878889, 9.852277768983186 ], [ 0.4661160492936929, 0.1380421246079127, 0.8738879793340638, -18.02322966460632 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 30, "182K", //(182_K_N, 182_K_CA, 182_K_C)
[ [ 0.49027649172895377, -0.8628024959452718, 0.12329158385128464, 158.8101446350395 ], [ 0.42434590709450826, 0.35986480572151547, 0.8309198954984498, -2.3151751526955824 ], [ -0.7612880616363874, -0.3550622122862769, 0.5425599622311683, 16.82757762893156 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 31, "183N", //(183_N_N, 183_N_CA, 183_N_C)
[ [ -0.4455452485592953, -0.6746232029558833, -0.5885345916084975, 175.7524548825811 ], [ -0.6141013263234852, -0.2480495253865145, 0.7492336044006978, 28.40321372030271 ], [ -0.6514361000850138, 0.6952373457957146, -0.30376971626045585, 31.909771239799547 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 32, "184W", //(184_W_N, 184_W_CA, 184_W_C)
[ [ 0.26417308132604084, 0.16703680292569933, -0.9499006735291149, 156.8575343328681 ], [ -0.8109801699286135, 0.5715576674250539, -0.1250319830690873, 51.70287137470065 ], [ 0.522038070532468, 0.8033806938656122, 0.28645368498021695, 8.47330870226984 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 33, "186T", //(186_T_N, 186_T_CA, 186_T_C)
[ [ 0.03810803541653977, -0.9983045716607932, -0.0439972708000397, 130.77325659555538 ], [ 0.2459254639266917, -0.0333056073772746, 0.968716368556696, 35.80054822302628 ], [ -0.9685393351997291, -0.047735926914548414, 0.2442393036603152, 53.80578233266481 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 34, "187E", //(187_E_N, 187_E_CA, 187_E_C)
[ [ -0.3715552741725895, -0.40654508064079414, -0.8346662660257105, 138.17011281471153 ], [ -0.8825800659062024, -0.124307575105236, 0.45343142153633287, 73.5583821027803 ], [ -0.28809565338535625, 0.9051346442260546, -0.3126214489158387, 52.567678810747026 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 35, "188T", //(188_T_N, 188_T_CA, 188_T_C)
[ [ 0.47733259202469297, -0.0851837585426929, -0.8745840862212897, 107.58847185493148 ], [ -0.6690432581224947, 0.6100193747768546, -0.4245674047282628, 81.1367327026563 ], [ 0.569679484755968, 0.7877944463357149, 0.23419050999299545, 31.236601785862987 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 36, "189L", //(189_L_N, 189_L_CA, 189_L_C)
[ [ 0.8852538967931776, -0.42078084274056426, -0.19816412539870742, 77.46838342813334 ], [ 0.3221636515439033, 0.8620443818185113, -0.39126725572048454, 58.0414709358094 ], [ 0.3354640365767881, 0.2825295845706801, 0.8986884410105241, 28.634639085622748 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 37, "190L", //(190_L_N, 190_L_CA, 190_L_C)
[ [ 0.08203838557208128, -0.952408900114266, -0.29357620863390915, 81.73355220749569 ], [ 0.6831310340844639, -0.160746457697933, 0.7123851252014305, 40.25732371868453 ], [ -0.7256732691531483, -0.258993944563769, 0.6374326969303361, 62.1171128590851 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 38, "191V", //(191_V_N, 191_V_CA, 191_V_C)
[ [ -0.07971221812049659, -0.5085621445286824, -0.8573275380125719, 80.19869571239911 ], [ -0.6726802207805076, -0.6072624441033063, 0.4227690203318667, 74.15894239927306 ], [ -0.7356271357508489, 0.6104071339147883, -0.2936934592614139, 79.34460516824416 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 39, "192Q", //(192_Q_N, 192_Q_CA, 192_Q_C)
[ [ 0.5005741371370113, -0.20070001686880676, -0.842107496973158, 54.11890185821118 ], [ -0.7236131970039669, 0.4369567951009724, -0.5342777370776225, 89.40717117263534 ], [ 0.4751941438519997, 0.876805715334901, 0.07350008982863578, 55.97120740195352 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 40, "193N", //(193_N_N, 193_N_CA, 193_N_C)
[ [ 0.7830600216818262, -0.5869339564060297, -0.2057317021297948, 30.86329151753662 ], [ 0.533273089723979, 0.8038444443538167, -0.2635411183435037, 59.39057209175528 ], [ 0.32005751704951213, 0.09665733339518774, 0.942454532421078, 59.089296789305514 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 41, "194A", //(194_A_N, 194_A_CA, 194_A_C)
[ [ -0.1341548770333732, -0.5776233900296591, -0.8052041283170404, 33.10291024191986 ], [ 0.3805853865559293, -0.7802667746816616, 0.4963250183780292, 58.13589103830502 ], [ -0.9149629678343459, -0.23986450262287964, 0.3245116143888292, 97.16465113195768 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 42, "195N", //(195_N_N, 195_N_CA, 195_N_C)
[ [ 0.9644912179933564, -0.15875276258697468, 0.21107877862233232, 2.7602781527535747 ], [ -0.1270320969152792, 0.42184590265681227, 0.8977242788100177, 67.72016566784843 ], [ -0.23155892720212506, -0.8926609629543923, 0.38669997213645885, 118.30735033745964 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 43, "196P", //(196_P_N, 196_P_CA, 196_P_C)
[ [ 0.17154352246682153, 0.7655775254231247, 0.6200515079142052, -2.576267140427336 ], [ -0.9841059687995644, 0.16249442760613242, 0.07163102100365641, 103.35133191837804 ], [ -0.04591581506132691, -0.6224842275624017, 0.7812842788404832, 132.20099754494697 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 44, "197D", //(197_D_N, 197_D_CA, 197_D_C)
[ [ -0.2090793420809778, 0.14964729226771764, 0.9663806272023078, 13.210670916336355 ], [ -0.07915524937996953, 0.9823890769256122, -0.16925172977798036, 96.0990880329633 ], [ -0.9746898353889715, -0.11188113985021762, -0.1935518931326577, 166.16588357861644 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 45, "198C", //(198_C_N, 198_C_CA, 198_C_C)
[ [ -0.8250684703350565, 0.3414396954309827, 0.4502010147067284, 44.414183530281164 ], [ 0.5589151604536219, 0.6101059512222995, 0.561591107210785, 78.22639933926118 ], [ -0.08292082171609425, 0.7149752881513675, -0.6942150060743406, 153.8158161704734 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 46, "199K", //(199_K_N, 199_K_CA, 199_K_C)
[ [ -0.2130190735932994, 0.8660620205921031, 0.4522824900141238, 49.25390804865761 ], [ -0.2378269858646122, -0.4949458148119359, 0.835743360843927, 99.47730210905371 ], [ 0.9476609093341384, 0.07046429512932326, 0.3114058188793549, 122.39829391318136 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 47, "200T", //(200_T_N, 200_T_CA, 200_T_C)
[ [ -0.09309140451833059, 0.4777852596271408, 0.8735303292318076, 55.94145878864247 ], [ -0.9678159681469466, 0.16260136969925817, -0.1920756266987398, 131.19835865367102 ], [ -0.233808131177239, -0.8632971911543271, 0.4472714137304733, 142.52957209335997 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 48, "201I", //(201_I_N, 201_I_CA, 201_I_C)
[ [ -0.585526637154647, 0.3201466432680621, 0.7447581379122772, 81.20140357471925 ], [ 0.2341826345157776, 0.9463447635778444, -0.22268830714710564, 113.40033162366016 ], [ -0.7760908779734627, 0.04401948719585599, -0.6290828513585422, 164.6916781491032 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 49, "202L", //(202_L_N, 202_L_CA, 202_L_C)
[ [ -0.5586977388137705, 0.7199172477217788, 0.4117960576269263, 99.22797870689087 ], [ 0.6534525456902788, 0.0763355549122439, 0.7531086598805871, 98.50885545546323 ], [ 0.5107412330870437, 0.6898492875178887, -0.5130802601307005, 134.62969457123052 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 50, "203K", //(203_K_N, 203_K_CA, 203_K_C)
[ [ -0.14258177800594254, 0.7629972839086079, 0.630480436912004, 102.6132062741406 ], [ -0.6513543946770914, -0.5519513691142548, 0.5206602910418074, 131.57141225058638 ], [ 0.7452569282572985, -0.33642953330673403, 0.5756798415816808, 116.06917850341708 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 51, "204A", //(204_A_N, 204_A_CA, 204_A_C)
[ [ -0.21646016785896954, 0.44075405346333285, 0.8711376814752696, 117.57138461353946 ], [ -0.7327885226485213, 0.5162559747863212, -0.4432840506628034, 146.6727789337134 ], [ -0.6451092750883871, -0.7343130346473953, 0.2112306567283, 147.75778270019168 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 52, "205L", //(205_L_N, 205_L_CA, 205_L_C)
[ [ -0.6536469055889931, 0.24693356746262687, 0.7153806931096995, 146.7217409404434 ], [ 0.6578786530914613, 0.6526480490203573, 0.3758273565302105, 122.00142137879267 ], [ -0.3740874237668229, 0.7162920754623086, -0.5890537004462325, 145.98349496944724 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 53, "206G", //(206_G_N, 206_G_CA, 206_G_C)
[ [ -0.5647727114358144, -0.7410227708824857, 0.3632038510961597, 163.14086350073265 ], [ -0.7042650597398069, 0.20337339307479094, -0.680183790323566, 139.49478556403642 ], [ 0.4301656774396635, -0.6399410254256948, -0.6367361886448624, 116.55908843500728 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 54, "207P", //(207_P_N, 207_P_CA, 207_P_C)
[ [ -0.1324945220283818, -0.7599308094622726, 0.636357106083122, 182.70842724028037 ], [ 0.3693400931375992, -0.6336325473471365, -0.6797776772912512, 124.1310263862093 ], [ 0.9198005747082246, 0.14496537439124546, 0.36462575744731174, 87.35177877426841 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 55, "208G", //(208_G_N, 208_G_CA, 208_G_C)
[ [ -0.8455729353028286, -0.4411383075243903, -0.3006715894757517, 202.35791522636353 ], [ 0.5294995160525233, -0.621169849343508, -0.5777354764654039, 91.82812842018181 ], [ 0.06809312434821978, -0.6477229437815832, 0.7588269331773997, 89.9896096847098 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 56, "209A", //(209_A_N, 209_A_CA, 209_A_C)
[ [ 0.33084608825881373, 0.9346337995090889, -0.13038606788704724, 180.55919173501528 ], [ 0.2700638000075076, -0.22616173573417825, -0.9359040619717433, 77.24784086294164 ], [ -0.904215908845529, 0.2744276409272728, -0.32723548109235595, 117.70775246507091 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 57, "210T", //(210_T_N, 210_T_CA, 210_T_C)
[ [ -0.28066233046818795, -0.32280652582879693, -0.903894132704991, 182.51502378196432 ], [ 0.8838368183130676, 0.2802830320387385, -0.3745315748310214, 39.291196688577784 ], [ 0.3742474246409708, -0.9040118189678239, 0.20664340375653684, 116.32500996759642 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 58, "211L", //(211_L_N, 211_L_CA, 211_L_C)
[ [ 0.8896999954765655, -0.4499909003298483, -0.07708506774552304, 152.9507786912195 ], [ 0.4261778407034761, 0.8791472611754584, -0.2132429160863562, 15.433637219640492 ], [ 0.16372649798466296, 0.15687027375517257, 0.9739534645298282, 115.38800702725761 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 59, "212E", //(212_E_N, 212_E_CA, 212_E_C)
[ [ -0.1334158066282124, -0.9562330869661917, -0.26041986470477113, 163.4590761998026 ], [ 0.5514150142869987, -0.28996302595642665, 0.7822166743281952, 5.086531973009405 ], [ -0.823493597158237, -0.03923935485326573, 0.5659669323114973, 150.56982971407805 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 60, "213E", //(213_E_N, 213_E_CA, 213_E_C)
[ [ -0.11977152410207288, -0.24327737767421972, -0.9625335835835661, 161.71281086475648 ], [ -0.803660749210649, -0.5454787290153208, 0.23787046132298284, 41.01655648654251 ], [ -0.5829100978645377, 0.8020405686146964, -0.13017966086823884, 163.3088459936197 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 61, "216T", //(216_T_N, 216_T_CA, 216_T_C)
[ [ -0.31144799421904357, -0.7499567614132088, -0.5835794743713542, 128.44844229127287 ], [ 0.014587536457042262, -0.6178249275529583, 0.7861803626867668, 23.79476091456424 ], [ -0.950151225162105, 0.23634131019522006, 0.2033603560632153, 196.45742235302467 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 62, "217A", //(217_A_N, 217_A_CA, 217_A_C)
[ [ 0.15375744607608702, 0.05151330695709095, -0.9867649299516595, 111.46842298975379 ], [ -0.9829489497482945, -0.0939472976417996, -0.1580672877433149, 57.57128710162158 ], [ -0.1008464672900543, 0.9942435740160193, 0.03618985442038097, 191.95081020179094 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 63, "218C", //(218_C_N, 218_C_CA, 218_C_C)
[ [ 0.9499731441580733, -0.1253581983407632, -0.2860705288686382, 77.06812051356664 ], [ -0.19762607797021783, 0.4679967973251364, -0.8613494824979723, 43.0657251189411 ], [ 0.24185731058732518, 0.874793872750703, 0.4198103399313707, 183.75712106600358 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 64, "219Q", //(219_Q_N, 219_Q_CA, 219_Q_C)
[ [ -0.09714447662552822, -0.5659566914107889, -0.8186916233286532, 77.63261943542044 ], [ 0.9659838626021235, -0.25169433090595333, 0.059372897706665616, 11.653536016220933 ], [ -0.23966252909754981, -0.7850751475095602, 0.5711557448799832, 205.4142195911702 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 65, "220G", //(220_G_N, 220_G_CA, 220_G_C)
[ [ 0.6104223151473642, 0.41137441693910765, -0.6768719866105712, 45.64840084762668 ], [ -0.27284532198897016, -0.6930466869025866, -0.6672643554409927, 1.1270827645273158 ], [ -0.7435993729414161, 0.5919944077956796, -0.31080957787651164, 223.0598327450092 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ]
   ]
 ]

];
