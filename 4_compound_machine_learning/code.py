from typing import List, Union
import numpy.typing as npt
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error


def draw_molecule(csvfile: str) -> None:
    # 課題 4-1
    df = pd.read_csv(csvfile)
    
    smiles = df[df['Compound ID'] == "CHEMBL540227"].iloc[0]['SMILES']
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol).save("data/4-1.png") # 保存先を指定するためにsaveメソッドを呼び出す
    
    pass

def create_2d_descriptors(smiles: str) -> Union[npt.NDArray[np.float_], List[float]]:
    # 課題 4-2
    mol = Chem.MolFromSmiles(smiles)
    desc_list = [desc_name[0] for desc_name in Descriptors._descList]
    desc_list.remove("AvgIpc")
    desc_list.remove("SPS")
    
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(desc_list)
    desc = calc.CalcDescriptors(mol)

    return desc

def predict_logpapp(csvfile: str) -> Union[npt.NDArray[np.float_], pd.Series, List[float]]:
    # 課題 4-3
    np.random.seed(0) # 出力を固定するためにseedを指定
    rfr = RandomForestRegressor(random_state=0) # 出力を固定するためにrandom_stateを指定

    df = pd.read_csv(csvfile)
    X = [create_2d_descriptors(smiles) for smiles in df['SMILES']]
    y = df['LogP app'].to_list()

    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=700, random_state=0)
    rfr.fit(X_train, y_train)
    pred = rfr.predict(X_test)

    return pred

def grid_search(csvfile: str) -> float:
    # 課題 4-4
    # こちらも出力を固定するためにseedやrandom_stateを指定すること
    np.random.seed(0)
    param = {
        "n_estimators": [100, 200, 400],
        "max_depth": [5,10,15]
    }
    rfr = RandomForestRegressor(random_state=0)

    df= pd.read_csv(csvfile)
    X = [create_2d_descriptors(smiles) for smiles in df['SMILES']]
    y = df['LogP app'].to_list()
    # Xは説明変数、yは目的変数
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=700, random_state=0)
    

    grid_search = GridSearchCV(rfr, param, cv=4, scoring="neg_mean_squared_error")
    grid_search.fit(X_train, y_train)
    y_pred = grid_search.predict(X_test)
    rmse = mean_squared_error(y_test, y_pred, squared=False)

    return rmse

if __name__ == "__main__":
    smiles = "C(=O)(c1ccc(OCCCCCC)cc1)CCNc1cc(Cl)ccc1"
    filepath = "data/fukunishi_data.csv"
    # 課題 4-1
    draw_molecule(filepath)
    # 課題 4-2
    print(create_2d_descriptors(smiles))
    # 課題 4-3
    print(predict_logpapp(filepath))
    # 課題 4-4
    print(grid_search(filepath))
