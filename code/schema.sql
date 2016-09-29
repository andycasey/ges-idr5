
/*
    Schema description for the iDR5 quality control phase.
*/

DROP TABLE IF EXISTS spectra;
CREATE TABLE spectra (
    cname char(16) not null,
    ges_fld char(12) not null,
    object char(32) not null,
    filename char(36) not null,
    ges_type char(8) not null,
    setup char(5) not null,
    wg integer not null,
    instrument char(18) not null,
    ra numeric not null,
    dec numeric not null,
    snr numeric not null,
    vel numeric,
    e_vel numeric,
    vrot numeric,
    e_vrot numeric,
    teff_irfm numeric,
    e_teff_irfm numeric,
    peculi char(11),
    remark char(11),
    tech char(69)
);
ALTER TABLE spectra ADD COLUMN id BIGSERIAL PRIMARY KEY;

DROP TABLE IF EXISTS recommended_idr4;
CREATE TABLE recommended_idr4 (
    cname char(16) not null,
    ges_fld char(23) not null,
    object char(74) not null,
    filename char(162) not null,
    ges_type char(20) not null,
    teff numeric,
    e_teff numeric,
    logg numeric,
    e_logg numeric,
    feh numeric,
    e_feh numeric,
    mh numeric,
    e_mh numeric,
    xi numeric,
    e_xi numeric,
    peculi char(29),
    remark char(16),
    tech char(87)
);
ALTER TABLE recommended_idr4 ADD COLUMN id BIGSERIAL PRIMARY KEY;

DROP TABLE IF EXISTS nodes;
CREATE TABLE nodes (
    wg integer not null,
    name char(10) not null
);
ALTER TABLE nodes ADD COLUMN id BIGSERIAL PRIMARY KEY;

DROP TABLE IF EXISTS results;
CREATE TABLE results (
    node_id integer not null,
    cname char(16) not null,
    filename char(255) not null,
    setup char(100) not null,
    
    snr numeric,

    vel numeric,
    e_vel numeric,
    vrot numeric,
    e_vrot numeric,
    
    teff numeric,
    e_teff numeric,
    nn_teff integer,
    enn_teff numeric,
    nne_teff numeric,
    sys_err_teff numeric,
    
    logg numeric,
    e_logg numeric,
    nn_logg integer,
    enn_logg numeric,
    nne_logg numeric,
    sys_err_logg numeric,
    lim_logg integer,
    
    feh numeric,
    e_feh numeric,
    nn_feh integer,
    enn_feh numeric,
    nne_feh numeric,
    sys_err_feh numeric,
    
    xi numeric,
    e_xi numeric,
    nn_xi integer,
    enn_xi numeric,
    nne_xi numeric,
    
    mh numeric,
    e_mh numeric,
    nn_mh integer,
    enn_mh numeric,
    nne_mh numeric,

    alpha_fe numeric,
    e_alpha_fe numeric,
    nn_alpha_fe integer,
    enn_alpha_fe numeric,
    nne_alpha_fe numeric,

    vrad numeric,
    e_vrad numeric,
    vsini numeric,
    e_vsini numeric,

    peculi char(255),
    remark char(255),
    tech char(255),

    propagated_peculi char(255),
    propagated_remark char(255),
    propagated_tech char(255),

    propagated_peculi_from_result_id integer,
    propagated_remark_from_result_id integer,
    propagated_tech_from_result_id integer,

    passed_quality_control boolean
);
ALTER TABLE results ADD COLUMN id BIGSERIAL PRIMARY KEY;
ALTER TABLE results ALTER COLUMN passed_quality_control SET DEFAULT true;


DROP TABLE IF EXISTS recommended_results;
CREATE TABLE recommended_results (
    wg integer not null,
    cname char(16) not null,
    provenance_setups char(255) not null,
    
    snr numeric,

    vel numeric,
    e_vel numeric,
    vrot numeric,
    e_vrot numeric,
    
    teff numeric,
    e_teff numeric,
    nn_teff integer,
    enn_teff numeric,
    nne_teff numeric,
    sys_err_teff numeric,
    
    logg numeric,
    e_logg numeric,
    nn_logg integer,
    enn_logg numeric,
    nne_logg numeric,
    sys_err_logg numeric,
    lim_logg integer,
    
    feh numeric,
    e_feh numeric,
    nn_feh integer,
    enn_feh numeric,
    nne_feh numeric,
    sys_err_feh numeric,
    
    xi numeric,
    e_xi numeric,
    nn_xi integer,
    enn_xi numeric,
    nne_xi numeric,
    
    mh numeric,
    e_mh numeric,
    nn_mh integer,
    enn_mh numeric,
    nne_mh numeric,

    alpha_fe numeric,
    e_alpha_fe numeric,
    nn_alpha_fe integer,
    enn_alpha_fe numeric,
    nne_alpha_fe numeric,

    vrad numeric,
    e_vrad numeric,
    vsini numeric,
    e_vsini numeric,

    peculi char(255),
    remark char(255),
    tech char(255)
);
ALTER TABLE recommended_results ADD COLUMN id BIGSERIAL PRIMARY KEY;
CREATE UNIQUE INDEX single_cname_result_per_wg ON recommended_results (wg, cname);

ALTER TABLE recommended_results ADD CONSTRAINT valid_e_teff_required
    CHECK ((e_teff > 0 AND e_teff is not null) OR teff = 'NaN');
ALTER TABLE recommended_results ADD CONSTRAINT valid_e_logg_required
    CHECK ((e_logg > 0 AND e_logg is not null) OR logg = 'NaN');
ALTER TABLE recommended_results ADD CONSTRAINT valid_e_feh_required
    CHECK ((e_feh > 0 AND e_feh is not null) OR feh = 'NaN');
ALTER TABLE recommended_results ADD CONSTRAINT valid_e_xi_required
    CHECK ((e_xi > 0 AND e_xi is not null) OR xi = 'NaN');
ALTER TABLE recommended_results ADD CONSTRAINT valid_e_mh_required
    CHECK ((e_mh > 0 AND e_mh is not null) OR mh = 'NaN');