
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
    teff numeric,
    e_teff numeric,
    logg numeric,
    e_logg numeric,
    mh numeric,
    e_mh numeric,
    xi numeric,
    e_xi numeric,
    peculi char(255),
    remark char(255),
    tech char(255)
);
ALTER TABLE results ADD COLUMN id BIGSERIAL PRIMARY KEY;
