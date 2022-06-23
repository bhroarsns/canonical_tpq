use std::fs::File;
use std::env;
use std::io::{BufReader, BufWriter, Read, Write};
use num::Complex;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    let args: Vec<String> = env::args().collect();

    let size = &args[1].trim().parse::<f64>().unwrap();
    let trial = *(&args[2].trim().parse::<usize>().unwrap());
    let l = &args[3].trim().parse::<f64>().unwrap();
    let here = &args[4];

    let temp_num = 80;
    let temp_delta = 0.05;

    let temps: Vec<f64> = (0..temp_num).map(|j| temp_delta * (temp_num - j) as f64).collect();
    let mut es: Vec<Vec<f64>> = (0..temp_num).map(|_| (0..trial).map(|_| 0.0).collect()).collect();
    let mut e2s: Vec<Vec<f64>> = (0..temp_num).map(|_| (0..trial).map(|_| 0.0).collect()).collect();

    for i in 0..trial {
        let output_name = &format!("{}gnuplot/ctpq_output_{}.dat", here, i);
        let mut output = BufWriter::new(File::create(output_name).expect("file not found: output"));

        let data: Vec<(Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>)> = {
            let ss_filename = &format!("{}output/SS_rand{}.dat", here, i);
            let norm_filename = &format!("{}output/Norm_rand{}.dat", here, i);

            let mut ss_f = BufReader::new(File::open(ss_filename).expect("file not found: ss"));
            let mut ss_data = String::new();
            ss_f.read_to_string(&mut ss_data).unwrap();

            let mut norm_f = BufReader::new(File::open(norm_filename).expect("file not found: norm"));
            let mut norm_data = String::new();
            norm_f.read_to_string(&mut norm_data).unwrap();
            
            ss_data.split('\n').zip(norm_data.split('\n')).flat_map(|(str1, str2)| {
                let dat1 = str1.split_whitespace().skip(1).take(2).map(|s| s.trim().parse::<f64>()).collect::<Result<Vec<f64>, _>>();
                match str2.split_whitespace().skip(1).next() {
                    Some(data) => {
                        match dat1 {
                            Ok(vec) => {
                                if vec.len() == 2 {
                                    let dat2 = data.trim().parse::<f64>();
                                    match dat2 {
                                        Ok(norm) => {
                                            let e = vec[0];
                                            let ee = vec[1];
        
                                            Ok((
                                                Complex::new(e / size, 0.0), //khk
                                                Complex::new(e * l / size - ee / size / size, 0.0), //khk1
                                                Complex::new(norm, 0.0), //kk
                                                Complex::new(l - e / size, 0.0) //kk1
                                            ))
                                        }
                                        Err(_) => Err(())
                                    }
                                } else {
                                    Err(())
                                }
                            }
                            Err(_) => Err(())
                        }
                    },
                    None => Err(())
                }
            }).collect()
        };

        for (ii, temp) in temps.iter().enumerate() {

            let (ene, nor): (Vec<[Complex<f64>;2]>, Vec<[Complex<f64>;2]>) = data.iter().enumerate().scan(
                Complex::new(0.0, 0.0),
                |factor, (k, (khk, khk1, kk, kk1))| {
                    let new_factor = *factor + if k == 0 {
                        Complex::new(0.0, 0.0)
                    } else {
                        2.0 * kk.ln()
                    };
                    let new_new_factor = new_factor + size.ln() - temp.ln() - ((2 * k + 1) as f64).ln();

                    *factor = new_new_factor + size.ln() - temp.ln() - ((2 * k + 2) as f64).ln();

                    Some(
                        (
                            [khk.ln() + new_factor, khk1.ln() + new_new_factor],
                            [new_factor, kk1.ln() + new_new_factor]
                        )
                    )
                }
            ).unzip();

            let mut ene: Vec<Complex<f64>> = ene.into_iter().flatten().collect();
            let mut nor: Vec<Complex<f64>> = nor.into_iter().flatten().collect();

            ene.sort_unstable_by(|a, b| b.re.partial_cmp(&a.re).unwrap());
            nor.sort_unstable_by(|a, b| b.re.partial_cmp(&a.re).unwrap());

            let divider = ene[0].re.min(nor[0].re);

            ene.sort_unstable_by(|a, b| a.re.partial_cmp(&b.re).unwrap());
            nor.sort_unstable_by(|a, b| a.re.partial_cmp(&b.re).unwrap());

            let ene = ene.into_iter().map(|x| {
                (x - divider).exp()
            }).sum::<Complex<f64>>();

            let nor = nor.into_iter().map(|x| {
                (x - divider).exp()
            }).sum::<Complex<f64>>();

            let ene = (ene/nor).re;

            writeln!(output, "{} {}", temp, ene).unwrap();

            es[ii][i] = ene;
            e2s[ii][i] = ene * ene;
        }
    };

    let varfilename = &format!("{}gnuplot/var.dat", here);
    let mut varfile = BufWriter::new(File::create(varfilename).unwrap());


    temps.into_iter().zip(
        es.into_iter().map(|x| {
            x.iter().sum::<f64>()
        }).zip(
            e2s.into_iter().map(|x| {
                x.iter().sum::<f64>()
            })
        )
    ).for_each(move |(t, (e, e2))| {
        writeln!(
            varfile,
            "{} {} {} {}",
            t,
            e / trial as f64,
            (e2 / (trial - 1) as f64 - e * e / trial as f64 / (trial - 1) as f64).sqrt(),
            e2 / (trial - 1) as f64 - e * e / trial as f64 / (trial - 1) as f64
        ).unwrap();
    });

    Ok(())
}