
/*
    Schema description for the iDR5 quality control phase.
*/

DROP TABLE IF EXISTS masterlist;
CREATE TABLE masterlist (
    cname char(21) not null,
    ges_fld char(23) not null,
    object char(28) not null,
    filename char(255) not null,
    ges_type char(8) not null,
    setup char(9) not null,
    wg integer not null,
    instrument char(1) not null,
    ra numeric not null,
    dec numeric not null,
    snr numeric not null,
    vel numeric,
    e_vel numeric,
    vrot numeric,
    e_vrot numeric,
    teff_irfm numeric,
    e_teff_irfm numeric,
    peculi char(255),
    remark char(255),
    tech char(255)
);
ALTER TABLE masterlist ADD COLUMN id BIGSERIAL PRIMARY KEY;

DROP TABLE IF EXISTS node;
CREATE TABLE node (
    wg integer not null,
    name char(10) not null
);
ALTER TABLE node ADD COLUMN id BIGSERIAL PRIMARY KEY;

DROP TABLE IF EXISTS result;
CREATE TABLE result (
    node_id integer not null,
    cname char(21) not null,
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
ALTER TABLE result ADD COLUMN id BIGSERIAL PRIMARY KEY;
