# Copyright (c) 2021 Hyejin Lee (KAIST)

import pandas as pd
import matplotlib.pyplot as plt

def main():
    colors = ["maroon", "tomato", "saddlebrown", "darkorange", "goldenrod", "darkolivegreen", "darkcyan", "dodgerblue", "rebeccapurple"]
    color_index = 0
    df = pd.read_csv("mutation_list.csv")
    df = df.transpose()
    df.columns = df.iloc[0]
    df = df.drop(df.index[0])
    df = df.reset_index()
    df['index'] = pd.to_datetime("2021-"+df['index'])

    plt.figure(figsize=(8,8))
    for mut in df.columns[1:]:
        if df[mut][5] < 50:
            plt.plot(df['index'], df[mut], color = "gray")
        else: 
            plt.plot(df['index'], df[mut], color = colors[color_index], label = mut)
            color_index = (color_index +1) % len(colors)
    plt.legend()
    plt.title("Occurence Plot")
    plt.savefig("task2_mutations.png")


    df2 = pd.read_csv("mutation_cumulative_list.csv")
    df2 = df2.transpose()
    df2.columns = df2.iloc[0]
    df2 = df2.drop(df2.index[0])
    df2 = df2.reset_index()
    df2['index'] = pd.to_datetime("2021-"+df2['index'])

    plt.figure(figsize=(8,8))
    for mut in df2.columns[1:]:
        if df2[mut][5] - df2[mut][0] < 20:
            plt.plot(df2['index'], df2[mut], color = "gray")
        else: 
            plt.plot(df2['index'], df2[mut], color = colors[color_index], label = mut)
            color_index = (color_index +1) % len(colors)
    plt.legend()
    plt.title("Cumulative Plot")
    plt.savefig("task2_mutations_cumulative.png")    

if __name__ == "__main__":
    main()

