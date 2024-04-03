#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("egg_nog", help="EggnogTsvFile")

    args = parser.parse_args()
    base_name = args.egg_nog.replace("-AMAPmastertable.tsv","")
    # Leer el archivo TSV
    df = pd.read_csv(args.egg_nog, sep='\t')
    df.replace("-", np.nan, inplace=True)
    # Suponiendo que la primera columna se llama 'primera_col', la dividimos
    df[['lineage', 'genome']] = df.iloc[:, 0].str.split('_', n=2, expand=True)[[0, 1]]

    total_conteo = df.shape[0]
    def contar_ocurrencias(row):
        if pd.notna(row['COG_category']):  # Verifica si el valor no es NaN
            for letra in row['COG_category']:
                yield row['lineage'], row['genome'], letra
        else:
            yield row['lineage'], row['genome'], 'NaN'

    ocurrencias = pd.DataFrame(
        [item for _, row in df.iterrows() for item in contar_ocurrencias(row)],
        columns=['lineage', 'genome', 'COG_category']
    )

    # Contamos las ocurrencias
    conteo = ocurrencias.groupby(['lineage', 'genome', 'COG_category']).size().reset_index(name='count')
    conteo_funcional = df.groupby(['lineage', 'genome', 'Functional_Category'], dropna=False).size().reset_index(name='count')    # Mostramos el resultado
    # Calcular el total de conteos por 'lineage' y 'genome' para el DataFrame conteo

    # Calcular el porcentaje
    conteo['percentage'] = (conteo['count'] / total_conteo) * 100

    # Repetir el proceso para el DataFrame conteo_funcional
    total_conteo_funcional = conteo_funcional.groupby(['lineage', 'genome'])['count'].transform('sum')
    conteo_funcional['percentage'] = (conteo_funcional['count'] / total_conteo_funcional) * 100

    # Mostrando el resultado (opcional)
    print(conteo.head())
    print(conteo_funcional.head())

    # Guardar los resultados en archivos TSV
    conteo.to_csv(base_name + '.AMAP.COG.tsv', sep='\t', index=False)
    conteo_funcional.to_csv(base_name + '.AMAP.MycoBrowser.tsv', sep='\t', index=False)



    # Supongamos que df1 y df2 son tus DataFrames
    # Aquí, reemplaza df1 y df2 con tus DataFrames reales y asegúrate de que las columnas se llamen 'categorias' y 'conteo'

    # Crear mapas de colores para las categorías en ambos DataFrames
    categorias_conteo = conteo['COG_category'].unique()
    colores_conteo = plt.cm.viridis(np.linspace(0, 1, len(categorias_conteo)))
    mapa_colores_conteo = dict(zip(categorias_conteo, colores_conteo))

    categorias_conteo_funcional = conteo_funcional['Functional_Category'].unique()
    colores_conteo_funcional = plt.cm.plasma(np.linspace(0, 1, len(categorias_conteo_funcional)))
    mapa_colores_conteo_funcional = dict(zip(categorias_conteo_funcional, colores_conteo_funcional))

    # Configurar el gráfico con dos subplots verticales
    fig, axs = plt.subplots(2, 1, figsize=(10, 12))  # 2 filas, 1 columna

    # Gráfico para el primer DataFrame (conteo)
    for categoria in categorias_conteo:
        axs[0].bar(categoria, conteo[conteo['COG_category'] == categoria]['percentage'], 
                   color=mapa_colores_conteo[categoria])
    axs[0].set_title('Eggnog mapper')
    axs[0].set_xlabel('')
    axs[0].set_ylabel('Count')
    axs[0].tick_params(axis='x', rotation=45)

    # Gráfico para el segundo DataFrame (conteo_funcional)
    for categoria in categorias_conteo_funcional:
        axs[1].bar(categoria, conteo_funcional[conteo_funcional['Functional_Category'] == categoria]['percentage'], 
                   color=mapa_colores_conteo_funcional[categoria])
    axs[1].set_title('MycoBrowser')
    axs[1].set_xlabel('')
    axs[1].set_ylabel('Count')
    axs[1].tick_params(axis='x', rotation=45)
    axs[1].set_xticklabels(axs[1].get_xticklabels(), rotation=45, fontsize=7, ha="right")
    plt.subplots_adjust(bottom=0.15, hspace=0.3)

    plt.savefig(base_name + '.AMAP-annot.png', dpi=300)

if __name__ == "__main__":
    main()