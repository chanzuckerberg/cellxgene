import React, { PureComponent } from "react";
import { connect, shallowEqual } from "react-redux";
import { Drawer, H3, H1, UL, Classes } from "@blueprintjs/core";
import Async from "react-async";
import {
  selectableCategoryNames,
  createCategorySummaryFromDfCol,
} from "../../util/stateManager/controlsHelpers";

@connect((state) => {
  return {
    annoMatrix: state.annoMatrix,
    schema: state.annoMatrix.schema,
    datasetTitle: state.config?.displayNames?.dataset ?? "",
    aboutURL: state.config?.links?.["about-dataset"],
    isOpen: state.controls.datasetDrawer,
    dataPortalProps: state.config?.["corpora_props"] ?? {},
  };
})
class InfoDrawer extends PureComponent {
  static watchAsync(props, prevProps) {
    return !shallowEqual(props.watchProps, prevProps.watchProps);
  }

  fetchAsyncProps = async (props) => {
    const { schema } = props.watchProps;
    const { annoMatrix } = this.props;

    const allCategoryNames = selectableCategoryNames(schema).sort();

    const nonUserAnnoCategories = allCategoryNames.map((catName) => {
      const isUserAnno = schema?.annotations?.obsByName[catName]?.writable;
      if (!isUserAnno) return annoMatrix.fetch("obs", catName);
      return null;
    });
    const singleValueCategories = (
      await Promise.all(nonUserAnnoCategories)
    ).reduce((acc, categoryData, i) => {
      const catName = allCategoryNames[i];

      const column = categoryData.icol(0);
      const colSchema = schema.annotations.obsByName[catName];

      const categorySummary = createCategorySummaryFromDfCol(column, colSchema);

      const { numCategoryValues } = categorySummary;
      //  Add to the array if the category has only one value
      if (numCategoryValues === 1) {
        acc.set(catName, categorySummary.allCategoryValues[0]);
      }
      return acc;
    }, new Map());

    return { singleValueCategories };
  };

  handleClose = () => {
    const { dispatch } = this.props;

    dispatch({ type: "toggle dataset drawer" });
  };

  render() {
    const {
      position,
      aboutURL,
      datasetTitle,
      schema,
      isOpen,
      dataPortalProps,
    } = this.props;

    return (
      <Drawer
        title="Dataset Overview"
        onClose={this.handleClose}
        {...{ isOpen, position }}
      >
        <Async
          watchFn={InfoDrawer.watchAsync}
          promiseFn={this.fetchAsyncProps}
          watchProps={{ schema }}
        >
          <Async.Pending>
            <InfoFormat
              skeleton
              {...{ datasetTitle, aboutURL, dataPortalProps }}
            />
          </Async.Pending>
          <Async.Rejected>
            {(error) => {
              console.error(error);
              return <span>Failed to load info</span>;
            }}
          </Async.Rejected>
          <Async.Fulfilled>
            {(asyncProps) => {
              const { singleValueCategories } = asyncProps;
              return (
                <InfoFormat
                  {...{
                    datasetTitle,
                    aboutURL,
                    singleValueCategories,
                    dataPortalProps,
                  }}
                />
              );
            }}
          </Async.Fulfilled>
        </Async>
      </Drawer>
    );
  }
}

const NUM_CATEGORIES = 8;

const singleValueCategoriesPlaceholder = Array.from(Array(NUM_CATEGORIES)).map(
  (_, index) => {
    return [index, index];
  }
);

const renderDOILink = (type, doi, skeleton) => {
  return (
    doi && (
      <>
        <H3 className={skeleton ? Classes.SKELETON : null}>{type}</H3>
        <p className={skeleton ? Classes.SKELETON : null}>
          <a href={doi} target="_blank" rel="noopener">
            {doi}
          </a>
        </p>
      </>
    )
  );
};

const buildAffiliations = (contributors = []) => {
  const affiliations = [];
  contributors.forEach((contributor) => {
    const { institution } = contributor;
    if (affiliations.indexOf(institution) === -1) {
      affiliations.push(institution);
    }
  });
  return affiliations;
};

const renderContributor = (contributor, affiliations) => {
  const { email, name, institution } = contributor;

  return (
    <span key={name}>
      {name}
      {email && `(${email})`}
      <sup>{affiliations.indexOf(institution) + 1}</sup>
    </span>
  );
};

const InfoFormat = ({
  datasetTitle,
  singleValueCategories = new Map(singleValueCategoriesPlaceholder),
  aboutURL = "thisisabouthtelengthofaurl",
  dataPortalProps = singleValueCategories,
  skeleton = false,
}) => {
  if (dataPortalProps.corpora_schema_version === "1.0.0") {
    dataPortalProps = {};
  }
  const {
    title,
    publication_doi: doi,
    preprint_doi: preprintDOI,
    organism,
    contributors,
    project_links: projectLinks,
  } = dataPortalProps;

  const affiliations = buildAffiliations(contributors);

  return (
    <div style={{ margin: 24, overflow: "auto" }}>
      <H1 className={skeleton ? Classes.SKELETON : null}>
        {title ?? datasetTitle}
      </H1>
      {contributors &&
      false /* Temp removed contributor section to avoid publishing PII */ && (
          <>
            <H3 className={skeleton ? Classes.SKELETON : null}>Contributors</H3>
            <p className={skeleton ? Classes.SKELETON : null}>
              {contributors.map((contributor) =>
                renderContributor(contributor, affiliations)
              )}
            </p>
            {affiliations.length > 0 && (
              <>
                <H3 className={skeleton ? Classes.SKELETON : null}>
                  Affiliations
                </H3>
                <UL>
                  {affiliations.map((item, index) => (
                    <div
                      id={`#afil${index}`}
                      key={item}
                      className={skeleton ? Classes.SKELETON : null}
                    >
                      <sup>{index + 1}</sup>
                      {"  "}
                      {item}
                    </div>
                  ))}
                </UL>
              </>
            )}
          </>
        )}
      {renderDOILink("DOI", doi, skeleton)}
      {renderDOILink("Preprint DOI", preprintDOI, skeleton)}
      {organism && (
        <>
          <H3 className={skeleton ? Classes.SKELETON : null}>Organism</H3>
          <p className={skeleton ? Classes.SKELETON : null}>{organism}</p>
        </>
      )}
      {singleValueCategories.size > 0 && (
        <>
          <H3 className={skeleton ? Classes.SKELETON : null}>
            Dataset Metadata
          </H3>
          <UL>
            {Array.from(singleValueCategories).map((pair) => {
              if (!pair[1] || pair[1] === "") return null;
              return (
                <li
                  className={skeleton ? Classes.SKELETON : null}
                  key={pair[0]}
                >{`${pair[0]}: ${pair[1]}`}</li>
              );
            })}
          </UL>
        </>
      )}
      {projectLinks ? (
        <>
          <H3 className={skeleton ? Classes.SKELETON : null}>Project Links</H3>
          <UL>
            {projectLinks.map((link) => {
              if (link.link_type === "SUMMARY") return null;
              return (
                <li
                  key={link.link_name}
                  className={skeleton ? Classes.SKELETON : null}
                >
                  <a href={link.link_url} target="_blank" rel="noopener">
                    {link.link_name}
                  </a>
                </li>
              );
            })}
          </UL>
        </>
      ) : (
        aboutURL && (
          <>
            <H3 className={skeleton ? Classes.SKELETON : null}>More Info</H3>
            <p>
              <a
                className={skeleton ? Classes.SKELETON : null}
                href={aboutURL}
                target="_blank"
                rel="noopener"
              >
                {aboutURL}
              </a>
            </p>
          </>
        )
      )}
    </div>
  );
};
export default InfoDrawer;
