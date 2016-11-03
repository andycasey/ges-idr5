
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


DROP TABLE IF EXISTS wg_recommended_results;
CREATE TABLE wg_recommended_results (
    provenance_ids_for_teff integer[],
    provenance_ids_for_logg integer[],
    provenance_ids_for_feh integer[],
    provenance_ids_for_mh integer[],
    provenance_ids_for_xi integer[],
    provenance_ids_for_alpha_fe integer[],    

    wg integer not null,
    cname char(16) not null,
    
    snr numeric,

    vel numeric,
    e_vel numeric,
    vrot numeric,
    e_vrot numeric,
    
    teff numeric,
    e_teff numeric,
    e_pos_teff numeric,
    e_neg_teff numeric,
    nn_nodes_teff integer,
    nn_spectra_teff integer,
    enn_teff numeric,
    nne_teff numeric,
    sys_err_teff numeric,
    
    logg numeric,
    e_logg numeric,
    e_pos_logg numeric,
    e_neg_logg numeric,
    nn_nodes_logg integer,
    nn_spectra_logg integer,
    enn_logg numeric,
    nne_logg numeric,
    sys_err_logg numeric,
    lim_logg integer,
    
    feh numeric,
    e_feh numeric,
    e_pos_feh numeric,
    e_neg_feh numeric,
    nn_nodes_feh integer,
    nn_spectra_feh integer,
    enn_feh numeric,
    nne_feh numeric,
    sys_err_feh numeric,
    
    xi numeric,
    e_xi numeric,
    e_pos_xi numeric,
    e_neg_xi numeric,
    nn_nodes_xi integer,
    nn_spectra_xi integer,
    enn_xi numeric,
    nne_xi numeric,
    
    mh numeric,
    e_mh numeric,
    e_pos_mh numeric,
    e_neg_mh numeric,
    nn_nodes_mh integer,
    nn_spectra_mh integer,
    enn_mh numeric,
    nne_mh numeric,

    alpha_fe numeric,
    e_alpha_fe numeric,
    e_pos_alpha_fe numeric,
    e_neg_alpha_fe numeric,
    nn_nodes_alpha_fe integer,
    nn_spectra_alpha_fe integer,
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
ALTER TABLE wg_recommended_results ADD COLUMN id BIGSERIAL PRIMARY KEY;
CREATE UNIQUE INDEX single_cname_result_per_wg ON wg_recommended_results (wg, cname);

ALTER TABLE wg_recommended_results ALTER COLUMN nn_nodes_teff SET DEFAULT 0;
ALTER TABLE wg_recommended_results ALTER COLUMN nn_nodes_logg SET DEFAULT 0;
ALTER TABLE wg_recommended_results ALTER COLUMN nn_nodes_feh SET DEFAULT 0;
ALTER TABLE wg_recommended_results ALTER COLUMN nn_nodes_xi SET DEFAULT 0;
ALTER TABLE wg_recommended_results ALTER COLUMN nn_nodes_mh SET DEFAULT 0;
ALTER TABLE wg_recommended_results ALTER COLUMN nn_nodes_alpha_fe SET DEFAULT 0;

ALTER TABLE wg_recommended_results ALTER COLUMN nn_spectra_teff SET DEFAULT 0;
ALTER TABLE wg_recommended_results ALTER COLUMN nn_spectra_logg SET DEFAULT 0;
ALTER TABLE wg_recommended_results ALTER COLUMN nn_spectra_feh SET DEFAULT 0;
ALTER TABLE wg_recommended_results ALTER COLUMN nn_spectra_xi SET DEFAULT 0;
ALTER TABLE wg_recommended_results ALTER COLUMN nn_spectra_mh SET DEFAULT 0;
ALTER TABLE wg_recommended_results ALTER COLUMN nn_spectra_alpha_fe SET DEFAULT 0;


ALTER TABLE wg_recommended_results ADD CONSTRAINT teff_provenance_required
    CHECK (teff = 'NaN' OR provenance_ids_for_teff is not null);
ALTER TABLE wg_recommended_results ADD CONSTRAINT logg_provenance_required
    CHECK (logg = 'NaN' OR provenance_ids_for_logg is not null);
ALTER TABLE wg_recommended_results ADD CONSTRAINT feh_provenance_required
    CHECK (feh = 'NaN' OR provenance_ids_for_feh is not null);
ALTER TABLE wg_recommended_results ADD CONSTRAINT mh_provenance_required
    CHECK (mh = 'NaN' OR provenance_ids_for_mh is not null);
ALTER TABLE wg_recommended_results ADD CONSTRAINT xi_provenance_required
    CHECK (xi = 'NaN' OR provenance_ids_for_xi is not null);
ALTER TABLE wg_recommended_results ADD CONSTRAINT alpha_fe_provenance_required
    CHECK (alpha_fe = 'NaN' OR provenance_ids_for_alpha_fe is not null);    

ALTER TABLE wg_recommended_results ADD CONSTRAINT valid_e_teff_required
    CHECK ((e_teff > 0 AND e_teff is not null) OR teff = 'NaN');
ALTER TABLE wg_recommended_results ADD CONSTRAINT valid_e_logg_required
    CHECK ((e_logg > 0 AND e_logg is not null) OR logg = 'NaN');
ALTER TABLE wg_recommended_results ADD CONSTRAINT valid_e_feh_required
    CHECK ((e_feh > 0 AND e_feh is not null) OR feh = 'NaN');
ALTER TABLE wg_recommended_results ADD CONSTRAINT valid_e_xi_required
    CHECK ((e_xi > 0 AND e_xi is not null) OR xi = 'NaN');
ALTER TABLE wg_recommended_results ADD CONSTRAINT valid_e_mh_required
    CHECK ((e_mh > 0 AND e_mh is not null) OR mh = 'NaN');


ALTER TABLE wg_recommended_results ADD COLUMN lim_vsini integer;
ALTER TABLE wg_recommended_results ADD COLUMN teff_phot numeric;
ALTER TABLE wg_recommended_results ADD COLUMN e_teff_phot numeric;
ALTER TABLE wg_recommended_results ADD COLUMN teff_irfm numeric;
ALTER TABLE wg_recommended_results ADD COLUMN e_teff_irfm numeric;
ALTER TABLE wg_recommended_results ADD COLUMN fbol_irfm numeric;

ALTER TABLE wg_recommended_results ADD COLUMN spt char(8);
ALTER TABLE wg_recommended_results ADD COLUMN veil numeric;
ALTER TABLE wg_recommended_results ADD COLUMN e_veil numeric;
ALTER TABLE wg_recommended_results ADD COLUMN ew_li numeric;
ALTER TABLE wg_recommended_results ADD COLUMN lim_ew_li integer;
ALTER TABLE wg_recommended_results ADD COLUMN e_ew_li numeric;
ALTER TABLE wg_recommended_results ADD COLUMN ewc_li numeric;
ALTER TABLE wg_recommended_results ADD COLUMN lim_ewc_li integer;
ALTER TABLE wg_recommended_results ADD COLUMN e_ewc_li numeric;
ALTER TABLE wg_recommended_results ADD COLUMN ew_ha_acc numeric;
ALTER TABLE wg_recommended_results ADD COLUMN e_ew_ha_acc numeric;
ALTER TABLE wg_recommended_results ADD COLUMN ha10 numeric;
ALTER TABLE wg_recommended_results ADD COLUMN e_ha10 numeric;
ALTER TABLE wg_recommended_results ADD COLUMN ew_ha_chr numeric;
ALTER TABLE wg_recommended_results ADD COLUMN e_ew_ha_chr numeric;
ALTER TABLE wg_recommended_results ADD COLUMN fha_chr numeric;
ALTER TABLE wg_recommended_results ADD COLUMN e_fha_chr numeric;
ALTER TABLE wg_recommended_results ADD COLUMN fwzi numeric;
ALTER TABLE wg_recommended_results ADD COLUMN e_fwzi numeric;
ALTER TABLE wg_recommended_results ADD COLUMN ew_hb_chr numeric;
ALTER TABLE wg_recommended_results ADD COLUMN e_ew_hb_chr numeric;
ALTER TABLE wg_recommended_results ADD COLUMN fhb_chr numeric;
ALTER TABLE wg_recommended_results ADD COLUMN e_fhb_chr numeric;
ALTER TABLE wg_recommended_results ADD COLUMN log_mdot_acc numeric;
ALTER TABLE wg_recommended_results ADD COLUMN e_log_mdot_acc numeric;
ALTER TABLE wg_recommended_results ADD COLUMN log_l_acc numeric;
ALTER TABLE wg_recommended_results ADD COLUMN e_log_l_acc numeric;
ALTER TABLE wg_recommended_results ADD COLUMN gamma numeric;
ALTER TABLE wg_recommended_results ADD COLUMN e_gamma numeric;
ALTER TABLE wg_recommended_results ADD COLUMN convol numeric;
ALTER TABLE wg_recommended_results ADD COLUMN e_convol numeric;
ALTER TABLE wg_recommended_results ADD COLUMN m_alpha numeric;
ALTER TABLE wg_recommended_results ADD COLUMN m_grid char(1);
ALTER TABLE wg_recommended_results ADD COLUMN m_broad numeric;
ALTER TABLE wg_recommended_results ADD COLUMN m_loops integer;
ALTER TABLE wg_recommended_results ADD COLUMN m_name  char(1);
