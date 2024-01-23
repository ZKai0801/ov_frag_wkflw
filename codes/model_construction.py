"""
10-fold CV model construction
"""

import pandas as pd
import os
import numpy as np
from sklearn.model_selection import RepeatedStratifiedKFold
import multiprocess



def model_per_split(df, fold_index, train_index, test_index):
    import numpy as np
    import pandas as pd
    import numpy as np
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.feature_selection import VarianceThreshold
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.ensemble import StackingClassifier
    from sklearn.ensemble import RandomForestClassifier
    from sklearn import svm
    import xgboost as xgb
    import pickle
    
    print(f"starting {fold_index} modelling...")
    
    X = np.array(df)
    y = np.array(df['group'])
    
    def make_base_classifier():
        glm = Pipeline([('scale', StandardScaler()),
                        ('selector', VarianceThreshold()),
                        ('lr', LogisticRegression(solver="newton-cg",
                                                random_state=123,
                                                max_iter=10000))])

        gbm = Pipeline([('scale', StandardScaler()),
                        ('selector', VarianceThreshold()),
                        ('gbm', GradientBoostingClassifier(random_state=123))])

        rf = Pipeline([('scale', StandardScaler()),
                       ('selector', VarianceThreshold()),
                       ('gbm', RandomForestClassifier(random_state=123))])

        xgboost = Pipeline([('scale', StandardScaler()),
                            ('selector', VarianceThreshold()),
                            ('xgb', xgb.XGBClassifier())])
        
        svm_model = Pipeline([('scale', StandardScaler()),
                              ('selector', VarianceThreshold()),
                              ('svm', svm.SVC(kernel = "linear", probability = True))])

        classifier = StackingClassifier(estimators=[("glm", glm),
                                                    ("gbm", gbm),
                                                    ("rf", rf),
                                                    ("xgboost", xgboost),
                                                    ("svm", svm_model)], 
                                        final_estimator=LogisticRegression())
        
        return classifier
        
    fsr_classifier = make_base_classifier()
    fsd_classifier = make_base_classifier()
    bpm_classifier = make_base_classifier()
    edm_classifier = make_base_classifier()
    cnv_classifier = make_base_classifier()
    
    
    fsr = [i for i,j in enumerate(df.columns.tolist()) if j.startswith("FSR")]
    fsd = [i for i,j in enumerate(df.columns.tolist()) if j.startswith("FSD")]
    cnv = [i for i,j in enumerate(df.columns.tolist()) if j.startswith("CNV")]
    edm = [i for i,j in enumerate(df.columns.tolist()) if j.startswith("EDM")]
    bpm = [i for i,j in enumerate(df.columns.tolist()) if j.startswith("BPM")]
    
    fsr_classifier.fit(X[train_index][:, fsr], y[train_index])
    fsd_classifier.fit(X[train_index][:, fsd], y[train_index])
    bpm_classifier.fit(X[train_index][:, bpm], y[train_index])
    edm_classifier.fit(X[train_index][:, edm], y[train_index])
    cnv_classifier.fit(X[train_index][:, cnv], y[train_index])

    proba_fsr = fsr_classifier.predict_proba(X[train_index][:, fsr])
    proba_fsd = fsd_classifier.predict_proba(X[train_index][:, fsd])
    proba_bpm = bpm_classifier.predict_proba(X[train_index][:, bpm])
    proba_edm = edm_classifier.predict_proba(X[train_index][:, edm])
    proba_cnv = cnv_classifier.predict_proba(X[train_index][:, cnv])
    
    cancer_scores_train = np.column_stack((proba_fsr[:, [1]], proba_fsd[:, [1]], 
                                           proba_bpm[:, [1]], proba_edm[:, [1]], 
                                           proba_cnv[:, [1]]))
    
    proba_fsr = fsr_classifier.predict_proba(X[test_index][:, fsr])
    proba_fsd = fsd_classifier.predict_proba(X[test_index][:, fsd])
    proba_bpm = bpm_classifier.predict_proba(X[test_index][:, bpm])
    proba_edm = edm_classifier.predict_proba(X[test_index][:, edm])
    proba_cnv = cnv_classifier.predict_proba(X[test_index][:, cnv])
    
    pred_fsr = fsr_classifier.predict(X[test_index][:, fsr])
    pred_fsd = fsd_classifier.predict(X[test_index][:, fsd])
    pred_bpm = bpm_classifier.predict(X[test_index][:, bpm])
    pred_edm = edm_classifier.predict(X[test_index][:, edm])
    pred_cnv = cnv_classifier.predict(X[test_index][:, cnv])
    
    cancer_scores_test = np.column_stack((proba_fsr[:, [1]], proba_fsd[:, [1]], 
                                          proba_bpm[:, [1]], proba_edm[:, [1]], 
                                          proba_cnv[:, [1]]))
    
    second_model = LogisticRegression()
    second_model.fit(cancer_scores_train, y[train_index])
    

    proba = second_model.predict_proba(cancer_scores_test)
    pred = second_model.predict(cancer_scores_test)
    
    prediction = df.loc[test_index, ['sampleID', 'group']]
    prediction['fold_index'] = fold_index
    
    prediction['proba_fsr'] = proba_fsr[:, 0]
    prediction['pred_fsr'] = pred_fsr
    prediction['proba_fsd'] = proba_fsd[:, 0]
    prediction['pred_fsd'] = pred_fsd
    prediction['proba_cnv'] = proba_cnv[:, 0]
    prediction['pred_cnv'] = pred_cnv
    prediction['proba_edm'] = proba_edm[:, 0]
    prediction['pred_edm'] = pred_edm
    prediction['proba_bpm'] = proba_bpm[:, 0]
    prediction['pred_bpm'] = pred_bpm
    
    prediction['cancer_scores'] = proba[:, 0]
    prediction['pred'] = pred
    
    prediction.to_csv(f"./{fold_index}_predictions.csv", index=None)
    
    pickle.dump(fsr_classifier, open(f"./{fold_index}_fsr_model.sav", 'wb'))
    pickle.dump(fsd_classifier, open(f"./{fold_index}_fsd_model.sav", 'wb'))
    pickle.dump(bpm_classifier, open(f"./{fold_index}_bpm_model.sav", 'wb'))
    pickle.dump(edm_classifier, open(f"./{fold_index}_edm_model.sav", 'wb'))
    pickle.dump(cnv_classifier, open(f"./{fold_index}_cnv_model.sav", 'wb'))
    pickle.dump(second_model, open(f"./{fold_index}_final_model.sav", "wb"))
    print(f"complete {fold_index} modelling")
    
    

if __name__ == "__main__":
    os.chdir("\\\\172.16.11.242\\test_data\\fragmentome\\results2\\models_test\\")

    df = pd.read_csv("./feature_matrix.csv")
    clinics = pd.read_excel("\\\\172.16.11.242\\test_data\\fragmentome\\manuscript\\version4\\Supplementary_Tables.xlsx", 
                            sheet_name="S1")
    df = df.merge(clinics[['sampleID', 'group']], how = "right", on = "sampleID")
    df['group'] = [1 if i == "OV" else 0 for i in df['group'].tolist()]
    df.reset_index(inplace=True, drop=True)
    X = np.array(df)
    y = np.array(df['group'])

    # -------------------------------------------------------
    print("starting modelling...")
    rskf = RepeatedStratifiedKFold(n_splits=10,
                                   n_repeats=10,
                                   random_state=123)
    pool = multiprocess.Pool(10)
    for fold_index, (train_index, test_index) in enumerate(rskf.split(X, y)):
        pool.apply_async(model_per_split, args=(df, fold_index, train_index, test_index))

    pool.close()
    pool.join()
