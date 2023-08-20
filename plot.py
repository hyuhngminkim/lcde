import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import json

if __name__=="__main__":
    # Configure argument parser
    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--dataset', type=str, required=True)
    parser.add_argument('-s', '--save', type=bool, default=False)
    parser.add_argument('-e', '--extension', type=str, default='pdf')
    parser.add_argument('-c', '--color', type=str, default='dodgerblue')
    parser.add_argument('-w', '--width', type=float, default=0.7)

    # Configure arguments
    args = parser.parse_args()
    dataset = args.dataset
    do_save = args.save
    extension = args.extension
    plot_color = args.color
    plot_linewidth = args.width

    # Read parameters from json file
    file_path = "./parameters/" + dataset + "_plot_parameters.json"
    file = open(file_path, 'r')
    data = json.load(file)
    parameters = data["parameters"]
    file.close()       
    plot_title = dataset.split(sep='_')[0]

    # Configure plot
    plt.figure(figsize=(8,6))
    plt.ylim([0, 1])
    plt.xlabel("x")
    plt.ylabel("Cumulative Distribution")
    plt.title(plot_title)

    # Main plot loop
    for idx, p in enumerate(parameters):
        slope = p["slope"]
        theta = p["theta"]
        addend = p["addend"]
        intercept = p["intercept"]
        for i, s in enumerate(slope):
            x = np.linspace(theta[i], theta[i + 1], 100)
            if s != 0:
                y = addend[i] + (1 if s >= 0 else -1) * np.exp(s * x + intercept[i])
            else:
                y = addend[i] + intercept[i] * x
            plt.plot(x, y, plot_color, linewidth=plot_linewidth)
        if idx < len(parameters) - 1:
            start_theta = theta[-1]
            end_theta = parameters[idx + 1]["theta"][0]
            x = np.linspace(start_theta, end_theta, 100)
            y_val = 0 * x + addend[-1] + (1 if slope[-1] >= 0 else -1) * np.exp(slope[-1] * start_theta + intercept[-1])
            plt.plot(x, y_val, plot_color, linewidth=plot_linewidth)

    if do_save:
        save_name = "./plot/" + dataset + "." + str.lower(extension)
        plt.savefig(save_name)
    else :
        plt.show()