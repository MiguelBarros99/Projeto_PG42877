from __future__ import division
from scipy import stats
import pandas as pd
import os
from sklearn.model_selection import train_test_split
import numpy as np
from sklearn.metrics import make_scorer
from sklearn.metrics import matthews_corrcoef
from propythia import shallow_ml
os.environ['TF_XLA_FLAGS'] = '--tf_xla_enable_xla_devices'
from collections import Counter
from sklearn.preprocessing import LabelEncoder
from keras.utils import np_utils
from propythia.deep_ml import DeepML
from keras.models import load_model


from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC

##### Data Load ##########
df = pd.read_csv('textfiles/Fulldata.csv')
df.dropna()

def turn_single_label(column, data):
    l = []
    for ec_list in data[column]:
        ec_l = set(ec_list)
        l.append(ec_l)
    data['ec_single_label']=l
    data = data.loc[data['ec_single_label'].apply(len)<2,:]
    return data

# get all the 1º level ec numbers complete
def get_ec_1_level(data, single_label=True):
    # get all until the 1 level (everything until the last dot
    l = []
    for ec_list in data['EC']:
        ec_list = str(ec_list)
        ec_1 = [x.strip()[0] for x in ec_list.split(';') ]
        # [^,]* = as many non-dot characters as possible,
        # . = a dot
        l.append(list(set(ec_1)))
    data['ec_number1']=l
    if single_label:
        data = turn_single_label('ec_number1', data)
    else:
        pass
    counts = Counter(x for xs in data['ec_number1'] for x in set(xs))
    counts.most_common()
    df = pd.DataFrame.from_dict(counts, orient='index').reset_index()
    return data, df

data, ec_counts = get_ec_1_level(data=df, single_label=True)

print(ec_counts)

ec_number = data['ec_number1']
flat_list = [item for sublist in ec_number for item in sublist] #the EC are in list of lists
encoder = LabelEncoder()
fps_y = encoder.fit_transform(flat_list) # fazer encoding das classes possiveis

vector=[]
from ast import literal_eval
for seq in range(df.shape[0]):
    new_vector = literal_eval(df['Vetor'][seq])
    vector.append(new_vector)

#treino
fps_x = np.array(vector)
fps_x = stats.zscore(fps_x)

print(fps_x.shape)
print(fps_y.shape)

colnames=[]
for  linha in range(1,20):
    for col in range(1,20):
        colnames.append(str(linha) + '_' + str(col))

###### ML Models ######
X_train, X_test, y_train, y_test = train_test_split(fps_x, fps_y, test_size=0.25, random_state = 42,stratify=fps_y)


###### Simple Cross Validation ######

print('RFCV5')
RFCVcc=shallow_ml.ShallowML(X_train, X_test, y_train, y_test, report_name='rfcv5Full', columns_names=colnames)
scoreRFCVcc = RFCVcc.cross_val_score_model(model_name='rf',model=None,score=make_scorer(matthews_corrcoef),cv=None,random_state=1,n_jobs=5)
print(scoreRFCVcc)

print('RFCV10')
RFCVdez=shallow_ml.ShallowML(X_train, X_test, y_train, y_test, report_name='rfcv10Full', columns_names=colnames)
scoreRFCVdez = RFCVdez.cross_val_score_model(model_name='rf',model=None,score=make_scorer(matthews_corrcoef),cv=10,random_state=1,n_jobs=5)
print(scoreRFCVdez)

print('ONEvsALL')
clf = OneVsRestClassifier(SVC(C= 32.0, kernel= 'rbf'), n_jobs=5).fit(X_train, y_train)
score = clf.score(X_test, y_test)
print(score)

print('SVMCV10')
SVMCVdez = shallow_ml.ShallowML(X_train, X_test, y_train, y_test, report_name='svm10Full', columns_names=colnames)
scoreSVMCVdez = SVMCVdez.cross_val_score_model(model_name='svm',model=None,score=make_scorer(matthews_corrcoef),cv=10,random_state=1,n_jobs=-1)
print(scoreSVMCVdez)

print('SVMCV5')
SVMCVcc = shallow_ml.ShallowML(X_train, X_test, y_train, y_test, report_name='svm5Full', columns_names=colnames)
scoreSVMCVcc = SVMCVcc.cross_val_score_model(model_name='svm',model=None,score=make_scorer(matthews_corrcoef),cv=5,random_state=1,n_jobs=-1)
print(scoreSVMCVcc)

##### Otimização de Hiperparametros

print('KNN')
KNN = shallow_ml.ShallowML(X_train, X_test, y_train, y_test, report_name='knnFull', columns_names=colnames)
best_knn_model = KNN.train_best_model(model_name='knn',model=None, scaler=None,score=make_scorer(matthews_corrcoef),
                         cv=5, optType='gridSearch', param_grid={'clf__n_neighbors': [7]},
                         n_jobs=4,random_state=1, n_iter=15, refit=True)

scores, scores_per_class, cm, cm2 = KNN.score_testset()
print(scores)
print(scores_per_class)
print(cm)

ROCknn = KNN.plot_roc_curve(classifier=best_knn_model, ylim=(0.0, 1.00), xlim=(0.0, 1.0),
                       title='ROC KNN Full',
                       path_save='ROC_PLOT/knnFull', show=False)

#################
print('RFN')
RF = shallow_ml.ShallowML(X_train, X_test, y_train, y_test, report_name='rfFull', columns_names=colnames)
best_rf_model = RF.train_best_model(model_name='rf',model=None, scaler=None,score=make_scorer(matthews_corrcoef),
                         cv=5, optType='gridSearch', param_grid={},
                         n_jobs=5,random_state=1, n_iter=15, refit=True)

scores, scores_per_class, cm, cm2 = RF.score_testset()
print(scores)
print(scores_per_class)
print(cm)

FTIrf = RF.features_importances_plot(classifier=best_rf_model , model_name='rf', top_features=30,
                             column_to_plot=None,
                             show=True, path_save='Feature/rfFull',
                             title=None,
                             kind='barh', figsize=(9, 7), color='r', edgecolor='black')
ROCrf = RF.plot_roc_curve(classifier=best_rf_model , ylim=(0.0, 1.00), xlim=(0.0, 1.0),
                       title='ROC RF Full',
                       path_save='ROC_PLOT/rfFull', show=False)

# ##############
print('SVM')
SVM = shallow_ml.ShallowML(X_train, X_test, y_train, y_test, report_name='svmFull', columns_names=colnames)
best_svm_model = SVM.train_best_model(model_name='svm',model=None, scaler=None,score=make_scorer(matthews_corrcoef),
                         cv=5, optType='gridSearch', param_grid={'clf__C': [32.0], 'clf__kernel': ['rbf']},
                         n_jobs=5,random_state=1, n_iter=15, refit=True, probability=True)
print(best_svm_model)
scores, scores_per_class, cm, cm2 = SVM.score_testset()
print(scores)
print(scores_per_class)
print(cm)

ROCsvm = SVM.plot_roc_curve(classifier=best_svm_model , ylim=(0.0, 1.00), xlim=(0.0, 1.0),
                       title='ROC SVM Full',
                       path_save='ROC_PLOT/svmFull', show=False)

ValCur = SVM.plot_validation_curve('clf__C', [22, 27, 37,42],
                              classifier=best_svm_model,
                              cv=5,
                              score=make_scorer(matthews_corrcoef), title="Validation Curve",
                              xlab="parameter range", ylab="MCC", n_jobs=1, show=False,
                              path_save='Validation/valcurFull')

###### DL Models ######

# confirmar que nao ha colunas so de zeros ou iguais

from sklearn.feature_selection import VarianceThreshold
sel = VarianceThreshold(0)
transf = sel.fit_transform(vector)

# try deep learning
# divide dataset in train test validation
x_train_1, x_test, y_train_1, y_test = train_test_split(fps_x, fps_y, test_size=0.20, random_state=42,
                                                        stratify=fps_y, shuffle=True)
x_train, x_dval, y_train, y_dval = train_test_split(x_train_1, y_train_1, test_size=0.25, random_state=42,
                                                    stratify=y_train_1, shuffle=True)

input_dim = fps_x.shape[1]
final_units= len(np.unique(fps_y))
dl = DeepML(x_train, y_train, x_test, y_test,
            number_classes=final_units, problem_type='multiclass',
            x_dval=x_dval, y_dval=y_dval,
            model=None,
            epochs=400, batch_size=512, callbacks=None,
            reduce_lr=True, early_stopping=True, checkpoint=True, tensorboard=False,
            early_stopping_patience=70, reduce_lr_patience=50, reduce_lr_factor=0.2, reduce_lr_min=0.00001,
            path='saved_models/Full', report_name='DNNFull', verbose=1,  validation_split=0.1, shuffle=True, class_weights=None)
#
# dnn simple  cv = None optType=None
dnn_simple = dl.run_dnn_simple(
    input_dim=input_dim,
    optimizer='Adam',
    hidden_layers=(256,256, 128,128),
    dropout_rate=(0.1,0.1,0.1,0.1),
    batchnormalization=(True,),
    l1=1e-5, l2=1e-4,
    final_dropout_value=0.3,
    initial_dropout_value=0.0,
    loss_fun=None, activation_fun=None,
    cv=None, optType=None, param_grid=None, n_iter_search=15, n_jobs=1,
    scoring=make_scorer(matthews_corrcoef))

score_simple = dl.model_simple_evaluate()
scores, report, cm, cm2 = dl.model_complete_evaluate()

print(scores)
print(report)
print(cm)

dl.roc_curve(classifier=dnn_simple, ylim=(0.0, 1.00), xlim=(0.0, 1.0),
                  title='ROC DNN Full',
                  path_save='ROC_PLOT/DnnFull.png', show=False, batch_size=None)
