# -*- coding: utf-8 -*-

from prog.fetch_data.sql_run import *
#Using database tables Patient and Exam

pats = s.query(Patient).all() #2448 patients
exs = s.query(Exam).all() #5599 exams; contains 2341 unique values of pt_id

high_risk = []
outcome = []
reg = []
for pat in pats:
    myexs = [ex for ex in exs if ex.pt_id==pat.pt_id]
    if sum([e.sty_indicator_high_risk_yn for e in myexs])>0:
        high_risk.append(pat)
    if sum([e.mri_cad_status_txt!=None and e.mri_cad_status_txt!='Unknown' for e in myexs ]):
        outcome.append(pat)
    if sum([e.previous_exam_reason_int=="Routine screening" for e in myexs])>0:
        reg.append(pat)
print "Pats with at least one high-risk",len(high_risk) # number of patients with one or more visits whose reason is "high risk"
print "Pats with at least one benign or malignant CAD outcome",len(outcome)
print "Pats with at least one 'previous exam reason' as 'routine screening'",len(reg)
"""
Patients with at least one high-risk 778
Patients with at least one benign or malignant CAD outcome 630
"""

exam_stats = {'yn': [],
'brca_1_yn': [],
'brca_2_yn': [],
'brca_1_or_2_yn': [],
'at_yn': [],
'other_gene_yn': [],
'prior_high_risk_marker_yn': [],
'prior_personal_can_hist_yn': [],
'hist_of_mantle_rad_yn': [],
'fam_hist_yn': []}
for ex in exs:
    for key in exam_stats:
        if getattr(ex,'sty_indicator_high_risk_'+key):
            exam_stats[key].append(ex)
for key in exam_stats:
    print key,len(exam_stats[key])
"""    
Number of exams in each high risk indicator category:
prior_personal_can_hist_yn 448
prior_high_risk_marker_yn 117
fam_hist_yn 739
brca_2_yn 374
yn 2746
brca_1_or_2_yn 245
brca_1_yn 434
hist_of_mantle_rad_yn 15
other_gene_yn 10
at_yn 11
"""
cads = [e.mri_cad_status_txt for e in exs]
cad_types = unique(cads)
for cad_t in cad_types:
    print cad_t,cads.count(cad_t)
"""
Number of exams with each CAD status type
None 1632
Benign by assumption 1328
Benign by pathology 183
Malignant 136
Unknown 2320
"""