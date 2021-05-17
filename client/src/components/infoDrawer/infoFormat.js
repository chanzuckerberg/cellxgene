import { H3, HTMLTable, Classes, Tooltip } from "@blueprintjs/core";
import React from "react";

const ONTOLOGY_KEY = "ontology_term_id";
const COLLECTION_LINK_ORDER_BY = [
  "DOI",
  "DATA_SOURCE",
  "RAW_DATA",
  "PROTOCOL",
  "LAB_WEBSITE",
  "OTHER",
];

const buildCollectionLinks = (links) => {
  /*
  sort links by custom sort order, create view-friendly model of link types.
   */
  const sortedLinks = [...links].sort(sortCollectionLinks);
  return sortedLinks.map((link) => {
    return {
      name: buildLinkName(link),
      type: transformLinkTypeToDisplay(link.type),
      url: link.url,
    };
  });
};

const buildDatasetMetadata = (singleValueCategories, corporaMetadata) => {
  /*
  transform Corpora metadata and single value categories into sort and render-friendly format.
  @returns [{key, value, tip}]
   */
  const metadata = [
    ...transformCorporaMetadata(corporaMetadata),
    ...transformSingleValueCategoriesMetadata(singleValueCategories),
  ];
  metadata.sort(sortDatasetMetadata);
  return metadata;
};

const sortCollectionLinks = (l0, l1) => {
  /*
  sort collection links by custom order.
  TODO(cc) revisit - improve readability here 
   */
  return (
    COLLECTION_LINK_ORDER_BY.indexOf(l1.type) -
    COLLECTION_LINK_ORDER_BY.indexOf(l0.type)
  );
};

const buildLinkName = (link) => {
  /*
  determine name to display for collection link.
  TODO(cc) error handling
   */
  if (link.name) {
    return link.name;
  }
  if (link.type === "DOI") {
    return new URL(link.url).pathname.substring(1);
  }
  return new URL(link.url).host;
};

const sortDatasetMetadata = (m0, m1) => {
  /*
  sort metadata key value pairs by key - alpha, ascending
   */
  if (m0.key < m1.key) {
    return -1;
  }
  if (m0.key > m1.key) {
    return 1;
  }
  return 0;
};

const transformCorporaMetadata = (corporaMetadata) => {
  /*
  build array of view model objects from given Corpora metadata object.
  @returns [{key, value}] 
   */
  return Object.entries(corporaMetadata)
    .filter(([, value]) => {
      return value;
    })
    .map(([key, value]) => {
      return {
        key,
        value,
      };
    });
};

const transformSingleValueCategoriesMetadata = (singleValueCategories) => {
  /*
  build array of view model objects from given single value categories map, ignoring ontology terms or metadata
  without values. add ontology terms as tooltips of their corresponding values.
  @returns [{key, value, tip}] where tip is an optional ontology term for the category
   */
  return Array.from(singleValueCategories.entries())
    .filter(([key, value]) => {
      if (key.indexOf(ONTOLOGY_KEY) >= 0) {
        // skip ontology terms
        return false;
      }
      // skip metadata without values
      return value;
    })
    .map(([key, value]) => {
      const viewModel = { key, value: String(value) };
      // add ontology term as tool tip if specified
      const tip = singleValueCategories.get(`${key}_${ONTOLOGY_KEY}`);
      if (tip) {
        viewModel.tip = tip;
      }
      return viewModel;
    });
};

const transformLinkTypeToDisplay = (type) => {
  /*
  convert link type from upper snake case to title case
  TODO(cc) revisit approach here, maybe create enum-type mapping to avoid string concat inside loop?
   */
  const tokens = type.split("_");
  return tokens
    .map((token) => {
      return token.charAt(0) + token.slice(1).toLowerCase();
    })
    .join(" ");
};

const renderCollectionLinks = (collection) => {
  /*
  render collection contact and links.
  TODO(cc) handle case where there is no contact and no links?
   */
  const links = buildCollectionLinks(collection.links);
  return (
    <>
      <p>Collection</p>
      <HTMLTable condensed>
        <tbody>
          <tr>
            <td>Contact</td>
            <td>{renderCollectionContactLink(collection.contact)}</td>
          </tr>
          {links.map(({ name, type, url }, i) => {
            return (
              <tr {...{ key: i }}>
                <td>{type}</td>
                <td>
                  <a href={url} target="_blank" rel="noopener">
                    {name || "link"}
                  </a>
                </td>
              </tr>
            );
          })}
        </tbody>
      </HTMLTable>
    </>
  );
};

const renderCollectionContactLink = (contact) => {
  /*
  display collection contact's name with a link to their associated email.
   */
  if (!contact || (!contact.name && !contact.email)) {
    return null;
  }
  const { name, email } = contact;
  if (email) {
    return <a href={`mailto:${email}`}>{name}</a>;
  }
  return name;
};

const renderDatasetMetadata = (singleValueCategories, corporaMetadata) => {
  /*
  render dataset metadata, mix of meta from Corpora and attributes found in categorical field.
   */
  if (
    singleValueCategories.size === 0 &&
    Object.entries(corporaMetadata).length === 0
  ) {
    return null;
  }
  const metadata = buildDatasetMetadata(singleValueCategories, corporaMetadata);
  return (
    <>
      <p>Dataset</p>
      <HTMLTable condensed>
        <tbody>
          {metadata.map(({ key, value, tip }) => {
            return (
              <tr {...{ key }}>
                <td>{key}</td>
                <td>
                  <Tooltip content={tip} minimal disabled={!tip}>
                    {value}
                  </Tooltip>
                </td>
              </tr>
            );
          })}
        </tbody>
      </HTMLTable>
    </>
  );
};

const InfoFormat = React.memo(
  ({ collection, singleValueCategories, dataPortalProps = {} }) => {
    if (
      ["1.0.0", "1.1.0"].indexOf(
        dataPortalProps.version?.corpora_schema_version
      ) === -1
    ) {
      dataPortalProps = {};
    }
    const { organism } = dataPortalProps;

    return (
      <div className={Classes.DIALOG_BODY}>
        <div className={Classes.DIALOG_BODY}>
          <H3>{collection.name}</H3>
          <p>{collection.description}</p>
          {renderCollectionLinks(collection)}
          {renderDatasetMetadata(singleValueCategories, { organism })}
        </div>
      </div>
    );
  }
);

export default InfoFormat;
