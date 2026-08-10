# Building and publishing RUFUS containers

RUFUS has **one** container definition — the root [`Dockerfile`](../Dockerfile). There is no
separate Singularity `.def`; SIFs are produced from the Docker image via `apptainer`. Builds,
tests, and publishing are automated by [`.github/workflows/build-publish.yml`](../.github/workflows/build-publish.yml).

## Submodule: modified jellyfish

RUFUS's custom jellyfish fork lives in its own repository
([stefinfection/modified-jellyfish](https://github.com/stefinfection/modified-jellyfish)) and is
pinned here as a git submodule at `src/modifiedJellyfish`. The container build compiles it from
that submodule. **Clone RUFUS with the submodule**, or the build (and any local `cmake`) will fail:
```bash
git clone --recursive https://github.com/stefinfection/RUFUS.git
# already cloned without --recursive:
git submodule update --init --recursive
```
To change the jellyfish source: edit/commit/tag the modified-jellyfish repo, then in RUFUS
`cd src/modifiedJellyfish && git checkout <new-tag> && cd ../.. && git add src/modifiedJellyfish`
and commit the updated pointer.

## Tag mapping (CI)

| Git event                | Docker Hub tag(s)                         | Stage |
|--------------------------|-------------------------------------------|-------|
| push to `docker` / `dev` | `stefinfection/rufus:dev`                 | dev   |
| push to `main`           | `stefinfection/rufus:stage`               | stage |
| push git tag `v*`        | `stefinfection/rufus:<version>` + `:latest` | prod |

Every build runs [`tests/smoke_test.sh`](../tests/smoke_test.sh) inside the freshly built image
and only pushes if it passes.

## Cutting a production release

1. Bump `RUFUS_VERSION` in [`resources/globals.txt`](../resources/globals.txt), commit, and merge
   to `main`. The release guard in CI fails the build if the git tag and `RUFUS_VERSION` disagree,
   so the bump must be committed *before* tagging and must be on the commit you tag.
2. Tag and push. Read the version out of `globals.txt` rather than typing it, so the tag and the
   guard can never disagree — this is the same extraction the workflow performs:
   ```bash
   VERSION=$(grep -E '^RUFUS_VERSION=' resources/globals.txt | cut -d'"' -f2)
   echo "tagging $VERSION"
   git tag "$VERSION"
   git push origin "$VERSION"
   ```
3. CI then automatically:
   - builds the image, smoke-tests it, and pushes `:<version>` + `:latest` to Docker Hub;
   - builds a SIF from that image with `apptainer`;
   - publishes it as a **new version of the existing Zenodo concept record** (shared concept DOI).

## Staging / using an image on HPC (CHPC)

Pull a published image straight into a SIF — no `sudo`, no manual `.def` build:
```bash
bash singularity/pull_staged_image.sh stage                                    # staging image
bash singularity/pull_staged_image.sh "$(grep -E '^RUFUS_VERSION=' resources/globals.txt \
                                          | cut -d'"' -f2)"                    # current prod release
```
This writes `rufus_<tag>.sif` into the zenodo_images dir and runs the smoke test against it.

## Required GitHub repo secrets

Set these on the `stefinfection/RUFUS` repository (Settings → Secrets and variables → Actions):

| Secret                     | Purpose                                            |
|----------------------------|----------------------------------------------------|
| `DOCKERHUB_USERNAME`       | Docker Hub login                                   |
| `DOCKERHUB_TOKEN`          | Docker Hub access token (push)                     |
| `ZENODO_TOKEN`             | Zenodo token (`deposit:write` + `deposit:actions`) |
| `ZENODO_CONCEPT_RECORD_ID` | Concept (all-versions) record id of the RUFUS Zenodo record |

## Known follow-up: give the Zenodo asset a stable name

**Do this after the first `v*` tag has published successfully — not before.**

`README.md` is the only place that still has to carry a literal version, and the reason is the
Zenodo download URL. Both halves of that URL move every release: the record id, because each
release is a new version record, and the filename, because
[`scripts/ci/zenodo_upload.sh`](../scripts/ci/zenodo_upload.sh) uploads
`$(basename "$SIF_PATH")`, i.e. `rufus_<version>.sif`. So the README needs editing every release,
which is exactly the drift [`scripts/ci/check_doc_versions.sh`](../scripts/ci/check_doc_versions.sh)
currently exists to police.

The durable fix removes the problem rather than guarding it:

1. Upload the SIF under a **stable name** (`rufus.sif`) in `zenodo_upload.sh` — either instead of,
   or in addition to, the versioned name. The image self-identifies its version anyway, via the
   OCI `org.opencontainers.image.version` label and the runtime banner.
2. Point the README at the **concept record**, which always resolves to the latest version, giving
   a permanent URL of the form
   `https://zenodo.org/records/<ZENODO_CONCEPT_RECORD_ID>/files/rufus.sif`.
3. Drop the README rules from `check_doc_versions.sh`; the guard becomes unnecessary.

Two things to confirm when doing this. The README currently references **two different** record ids
(`18284901` in the prose link, and the download URL) — establish which is the concept record and use
that one consistently. And the reason for waiting is that step 1 modifies the release job, which is
the single irreversible, least-exercised part of the pipeline: a published Zenodo version cannot be
retracted, so it should not be the change under test on the run that first proves the job works.
