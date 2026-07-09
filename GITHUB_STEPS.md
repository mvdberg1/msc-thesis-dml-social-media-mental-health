# GitHub upload steps

The local thesis folder is already a Git repository on branch `main`, with a
first commit named `Prepare thesis replication package`.

## Option A: GitHub website plus Terminal

1. Go to <https://github.com> and sign in.
2. Click the `+` button in the top-right corner and choose `New repository`.
3. Use a clear repository name, for example:
   `msc-thesis-dml-social-media-mental-health`.
4. Add a short description:
   `Replication package for an MSc thesis on social media use, adolescent mental well-being, panel models, and Double Machine Learning.`
5. Choose `Public` if the Canvas marker must be able to open the link without
   being invited. Choose `Private` only if you will invite the marker or provide
   access separately.
6. Do not tick `Add a README file`.
7. Do not add a `.gitignore`.
8. Do not add a licence unless your supervisor explicitly asks for one.
9. Click `Create repository`.
10. Copy the HTTPS repository URL. It will look like:
    `https://github.com/YOUR-USERNAME/msc-thesis-dml-social-media-mental-health.git`.
11. In Terminal, run:

```sh
cd "/Users/m.vandenberg/Documents/MSc Econometrics/Thesis"
git remote add origin https://github.com/YOUR-USERNAME/msc-thesis-dml-social-media-mental-health.git
git push -u origin main
```

After the push, the repository page should show `README.md`, `thesis.pdf`,
`thesis.tex`, `R/`, `tables/`, `figures/`, and the replication documentation.

## Option B: GitHub Desktop

1. Open GitHub Desktop.
2. Choose `File` > `Add Local Repository...`.
3. Select:
   `/Users/m.vandenberg/Documents/MSc Econometrics/Thesis`.
4. Click `Add Repository`.
5. Click `Publish repository`.
6. Use a clear repository name, for example:
   `msc-thesis-dml-social-media-mental-health`.
7. Choose public/private visibility as described above.
8. Publish.

## After the repository is online

Copy the repository URL, for example:

```text
https://github.com/YOUR-USERNAME/msc-thesis-dml-social-media-mental-health
```

This URL can then be added to the thesis in the reproducibility paragraph at
the end of Chapter 5.
