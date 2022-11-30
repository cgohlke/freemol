struct t_field {
        char    name[11];
        int     type;
        float   bondunit, bond_cubic, bond_quartic;
        float   angleunit, angle_cubic, angle_quartic, angle_pentic, angle_sextic;
        float   str_bndunit, ang_angunit, str_torunit, torsionunit;
        char    vdwtype[14], radiusrule[11], radiustype[11], radiussize[11], epsrule[11];
        float   vdw_14scale, a_expterm, b_expterm, c_expterm;
        float   chg_14scale, dielectric;
        } field;

