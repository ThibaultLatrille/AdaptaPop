import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import scipy.stats as st
from cds_libraries import params_pb, load_pb_table, str_to_table, columns


def p_value(v, mean, std):
    if mean > v:
        return st.norm.cdf((v - mean) / std)
    else:
        return st.norm.cdf((mean - v) / std)

data_path = "/mnt/sda1/AdaptaPop/data"
pb_table = load_pb_table(data_path)

for negative in [False, True]:
    my_dpi = 96
    plt.figure(figsize=(2 * 1920 / my_dpi, 2 * 1440 / my_dpi), dpi=my_dpi)
    n_plot = 1
    txt_file = open('{0}/79_GRCh38_correlation_pb{1}.out'.format(data_path, "" if negative else "_only_pos"), 'w')
    for x, x_label in params_pb:
        for y, y_label in params_pb:
            plt.subplot(len(params_pb), len(params_pb), n_plot)
            n_plot += 1
            table_x = str_to_table(pb_table, x)
            table_y = str_to_table(pb_table, y)
            mask = [True if negative or (val_x > 0 and val_y > 0) else False for val_x, val_y in zip(table_x, table_y)]
            x_val = table_x[mask]
            y_val = table_y[mask]
            plt.scatter(x_val, y_val, linewidth=3, label=r"$n={0}$".format(len(x_val)))
            idf = np.linspace(min((min(x_val), 0)), max(x_val), 30)
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            plt.xlim((min((min(x_val), 0)), max(x_val)))
            plt.ylim((min((min(y_val), 0)), max(y_val)))
            model = sm.OLS(y_val, sm.add_constant(x_val))
            results = model.fit()
            b, a = results.params[0:2]
            plt.plot(idf, a * idf + b, 'r-',
                     label=r"$y={0:.3g}x + {1:.3g}$ ($r^2={2:.3g})$".format(float(a), float(b), results.rsquared))
            if x != y:
                prstd, iv_l, iv_u = wls_prediction_std(results)
                sorted_x = sorted(x_val)
                plt.plot(sorted_x, sorted(iv_u), 'r--')
                plt.plot(sorted_x, sorted(iv_l), 'r--')
                tr_ids = [i for i, y in enumerate(x_val) if (y < iv_l[i] or y > iv_u[i])]
                txt_file.write("x={0}; y={1}".format(x, y))
                txt_file.write("\nCdsId\tx_{0}\ty_{1}\ty_predicted\ty_lower_bound\ty_upper_bound\tp_value\t".format(x, y) +
                               "\t".join(columns) + "\n")
                txt_file.write("\n".join(["\t".join(
                    [str(j) for j in ([pb_table["CdsId"][mask][i][2:-1],
                                       x_val[i],
                                       y_val[i],
                                       results.fittedvalues[i],
                                       iv_l[i],
                                       iv_u[i],
                                       p_value(y_val[i], results.fittedvalues[i], prstd[i])] +
                                      [pb_table[col][mask][i] for col in columns])
                     ]) for i in tr_ids]) + "\n\n")
            plt.legend()

    txt_file.close()
    plt.tight_layout()
    plt.savefig("{0}/79_GRCh38_correlation_pb{1}.svg".format(data_path, "" if negative else "_only_pos"), format="svg")
    plt.savefig("{0}/79_GRCh38_correlation_pb{1}.png".format(data_path, "" if negative else "_only_pos"), format="png")

print('Plot completed')
