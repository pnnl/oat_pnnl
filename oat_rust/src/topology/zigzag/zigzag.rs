







pub fn zigzag() {
    let Vi;
    let Vii;
    let Ai;
    let Si;
    let C0 = identity;
    let j = i-1;
    let k = i+1;
    Z = [];

    for i in 0 .. N {
        // left basis
        // ----------
        if Si < 0 {
            Di = Ai * Cj;
            Ri * Mi = Di * Ci;
            Wi = Ri; 
        } else { 
            .. 
        }
        // right basis 
        // -----------       
        if Sk < 0 {
            Dk = Ak;
            Rk * Mk = Dk * Ck;
            Wk = Rk; 
        } else { 
            .. 
        }
        // zip bases
        // ---------
        Dl = Wi^{-1} * Wk;
        Rl * Ml = Dl * Cl;
        // store zipped basis
        // ------------------
        Z.append(Dl)
    }
    return Z
}