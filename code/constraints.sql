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