import { H3, H1, UL, Classes } from "@blueprintjs/core";
import React from "react";

import Truncate from "../util/truncate";

const renderContributors = (contributors, affiliations, skeleton) => {
  // eslint-disable-next-line no-constant-condition --  Temp removed contributor section to avoid publishing PII
  if (!contributors || contributors.length === 0 || true) return null;
  return (
    <>
      <H3 className={skeleton ? Classes.SKELETON : null}>Contributors</H3>
      <p className={skeleton ? Classes.SKELETON : null}>
        {contributors.map((contributor) => {
          const { email, name, institution } = contributor;

          return (
            <span key={name}>
              {name}
              {email && `(${email})`}
              <sup>{affiliations.indexOf(institution) + 1}</sup>
            </span>
          );
        })}
      </p>
      {renderAffiliations(affiliations, skeleton)}
    </>
  );
};

// generates a list of unique institutions by order of appearance in contributors
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

const renderAffiliations = (affiliations, skeleton) => {
  if (affiliations.length === 0) return null;
  return (
    <>
      <H3 className={skeleton ? Classes.SKELETON : null}>Affiliations</H3>
      <UL>
        {affiliations.map((item, index) => (
          <div key={item} className={skeleton ? Classes.SKELETON : null}>
            <sup>{index + 1}</sup>
            {"  "}
            {item}
          </div>
        ))}
      </UL>
    </>
  );
};

const renderDOILink = (type, doi, skeleton) => {
  if (!doi) return null;
  return (
    <>
      <H3 className={skeleton ? Classes.SKELETON : null}>{type}</H3>
      <p className={skeleton ? Classes.SKELETON : null}>
        <a href={doi} target="_blank" rel="noopener">
          {doi}
        </a>
      </p>
    </>
  );
};

const renderOrganism = (organism, skeleton) => {
  if (!organism) return null;
  return (
    <>
      <H3 className={skeleton ? Classes.SKELETON : null}>Organism</H3>
      <p className={skeleton ? Classes.SKELETON : null}>{organism}</p>
    </>
  );
};

const ONTOLOGY_TERM_ID = "ontology_term_id";
const CAT_WIDTH = "30%";
const VAL_WIDTH = "35%";
// Render list of metadata attributes found in categorical field
const renderSingleValueCategories = (singleValueCategories, skeleton) => {
  if (singleValueCategories.size === 0) return null;
  return (
    <>
      <H3 className={skeleton ? Classes.SKELETON : null}>Dataset Metadata</H3>
      <UL>
        {Array.from(singleValueCategories).map((pair, index, array) => {
          const [category, value] = pair;
          let ontologyTermID = "";
          // If the value is empty skip it
          if (!value || value === "") return null;

          // If the next category is a ontology term, let's save it's value to render now
          if (array[index + 1]?.[0]?.toString().includes(ONTOLOGY_TERM_ID))
            ontologyTermID = array[index + 1][1];
          // If this category is an ontology term skip it since we already used it
          if (category.toString().includes(ONTOLOGY_TERM_ID)) {
            return null;
          }

          // Create the list item, appending an ontology term if we found it
          return (
            <li
              className={skeleton ? Classes.SKELETON : null}
              key={category}
              style={{ width: "100%" }}
            >
              <Truncate>
                <span style={{ width: CAT_WIDTH }}>{`${category}:`}</span>
              </Truncate>
              <Truncate>
                <span style={{ width: VAL_WIDTH }}>{value}</span>
              </Truncate>
              <Truncate>
                <span style={{ width: VAL_WIDTH }}>{ontologyTermID}</span>
              </Truncate>
            </li>
          );
        })}
      </UL>
    </>
  );
};

// Renders any links found in the config where link_type is not "SUMMARY"
// If there are no links in the config, render the aboutURL
const renderLinks = (projectLinks, aboutURL, skeleton) => {
  if (!projectLinks && !aboutURL) return null;
  if (projectLinks)
    return (
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
    );

  return (
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
  );
};

const NUM_CATEGORIES = 8;

// Generates arbitrary placeholder array for singleValueCategories skeleton shape
const singleValueCategoriesPlaceholder = Array.from(Array(NUM_CATEGORIES)).map(
  (_, index) => {
    return [index, index];
  }
);

const InfoFormat = React.memo(
  ({
    datasetTitle,
    singleValueCategories = new Map(singleValueCategoriesPlaceholder),
    aboutURL = "thisisabouthtelengthofaurl",
    dataPortalProps = {},
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
        {renderContributors(contributors, affiliations, skeleton)}
        {renderDOILink("DOI", doi, skeleton)}
        {renderDOILink("Preprint DOI", preprintDOI, skeleton)}
        {renderOrganism(organism, skeleton)}
        {renderSingleValueCategories(singleValueCategories, skeleton)}
        {renderLinks(projectLinks, aboutURL, skeleton)}
      </div>
    );
  }
);

export default InfoFormat;
