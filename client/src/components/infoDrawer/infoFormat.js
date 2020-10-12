import { H3, H1, UL } from "@blueprintjs/core";
import React from "react";

import Truncate from "../util/truncate";

const renderContributors = (contributors, affiliations) => {
  // eslint-disable-next-line no-constant-condition --  Temp removed contributor section to avoid publishing PII
  if (!contributors || contributors.length === 0 || true) return null;
  return (
    <>
      <H3>Contributors</H3>
      <p>
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
      {renderAffiliations(affiliations)}
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

const renderAffiliations = (affiliations) => {
  if (affiliations.length === 0) return null;
  return (
    <>
      <H3>Affiliations</H3>
      <UL>
        {affiliations.map((item, index) => (
          <div key={item}>
            <sup>{index + 1}</sup>
            {"  "}
            {item}
          </div>
        ))}
      </UL>
    </>
  );
};

const renderDOILink = (type, doi) => {
  if (!doi) return null;
  return (
    <>
      <H3>{type}</H3>
      <p>
        <a href={doi} target="_blank" rel="noopener">
          {doi}
        </a>
      </p>
    </>
  );
};

const renderOrganism = (organism) => {
  if (!organism) return null;
  return (
    <>
      <H3>Organism</H3>
      <p>{organism}</p>
    </>
  );
};

const ONTOLOGY_KEY = "ontology_term_id";
const CAT_WIDTH = "30%";
const VAL_WIDTH = "35%";
// Render list of metadata attributes found in categorical field
const renderSingleValueCategories = (singleValueCategories) => {
  if (singleValueCategories.size === 0) return null;
  return (
    <>
      <H3>Dataset Metadata</H3>
      <UL>
        {Array.from(singleValueCategories).reduce((elems, pair) => {
          const [category, value] = pair;
          // If the value is empty skip it
          if (!value) return elems;

          // If this category is a ontology term, let's add its value to the previous node
          if (String(category).includes(ONTOLOGY_KEY)) {
            const prevElem = elems.pop();
            // Props aren't extensible so we must clone and alter the component to append the new child
            elems.push(
              React.cloneElement(
                prevElem,
                prevElem.props,
                // Concat returns a new array
                prevElem.props.children.concat([
                  <Truncate key="ontology">
                    <span style={{ width: VAL_WIDTH }}>{value}</span>
                  </Truncate>,
                ])
              )
            );
          } else {
            // Create the list item
            elems.push(
              <li key={category} style={{ width: "100%" }}>
                <Truncate>
                  <span style={{ width: CAT_WIDTH }}>{`${category}:`}</span>
                </Truncate>
                <Truncate>
                  <span style={{ width: VAL_WIDTH }}>{value}</span>
                </Truncate>
              </li>
            );
          }
          return elems;
        }, [])}
      </UL>
    </>
  );
};

// Renders any links found in the config where link_type is not "SUMMARY"
// If there are no links in the config, render the aboutURL
const renderLinks = (projectLinks, aboutURL) => {
  if (!projectLinks && !aboutURL) return null;
  if (projectLinks)
    return (
      <>
        <H3>Project Links</H3>
        <UL>
          {projectLinks.map((link) => {
            if (link.link_type === "SUMMARY") return null;
            return (
              <li key={link.link_name}>
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
      <H3>More Info</H3>
      <p>
        <a href={aboutURL} target="_blank" rel="noopener">
          {aboutURL}
        </a>
      </p>
    </>
  );
};

const InfoFormat = React.memo(
  ({
    datasetTitle,
    singleValueCategories,
    aboutURL = "thisisabouthtelengthofaurl",
    dataPortalProps = {},
  }) => {
    if (dataPortalProps.version?.corpora_schema_version !== "1.0.0") {
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
        <H1>{title ?? datasetTitle}</H1>
        {renderContributors(contributors, affiliations)}
        {renderDOILink("DOI", doi)}
        {renderDOILink("Preprint DOI", preprintDOI)}
        {renderOrganism(organism)}
        {renderSingleValueCategories(singleValueCategories)}
        {renderLinks(projectLinks, aboutURL)}
      </div>
    );
  }
);

export default InfoFormat;
