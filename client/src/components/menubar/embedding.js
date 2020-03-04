import React from "react";
import {
	AnchorButton,
	ButtonGroup,
	Popover,
	Button,
	Radio,
	RadioGroup,
	Tooltip,
	Position
} from "@blueprintjs/core";
import { connect } from "react-redux";
import * as globals from "../../globals";
import { World } from "../../util/stateManager";
import actions from "../../actions";

@connect(state => ({
	universe: state.universe,
	world: state.world,
	layoutChoice: state.layoutChoice,
	reembedController: state.reembedController,
	enableReembedding: state.config?.parameters?.["enable-reembedding"] ?? false
}))
class Embedding extends React.PureComponent {
	handleLayoutChoiceChange = e => {
		const { dispatch } = this.props;
		dispatch({
			type: "set layout choice",
			layoutChoice: e.currentTarget.value
		});
	};

	renderReembedding() {
		const {
			enableReembedding,
			world,
			universe,
			dispatch,
			reembedController
		} = this.props;

		if (!enableReembedding) return null;

		const loading = !!reembedController?.pendingFetch;
		const disabled = World.worldEqUniverse(world, universe);
		const tipContent = disabled
			? "Subset cells first, then click to recompute UMAP embedding."
			: "Click to recompute UMAP embedding on the current cell subset.";

		return (
			<Tooltip
				content={tipContent}
				position="bottom"
				hoverOpenDelay={globals.tooltipHoverOpenDelay}
			>
				<AnchorButton
					icon="new-object"
					style={{ marginRight: 10 }}
					disabled={disabled}
					onClick={() => dispatch(actions.requestReembed())}
					loading={loading}
				/>
			</Tooltip>
		);
	}

	render() {
		const { layoutChoice } = this.props;

		return (
			<ButtonGroup
				style={{
					marginRight: 10
				}}
			>
				<Popover
					target={
						<Tooltip
							content="Select embedding for visualization"
							position="bottom"
							hoverOpenDelay={globals.tooltipHoverOpenDelay}
						>
							<Button
								type="button"
								data-testid="layout-choice"
								icon="heatmap"
								style={{
									cursor: "pointer"
								}}
							/>
						</Tooltip>
					}
					position={Position.BOTTOM_RIGHT}
					content={
						<div
							style={{
								display: "flex",
								justifyContent: "flex-start",
								alignItems: "flex-start",
								flexDirection: "column",
								padding: 10
							}}
						>
							<RadioGroup
								label="Embedding Choice"
								onChange={this.handleLayoutChoiceChange}
								selectedValue={layoutChoice.current}
							>
								{layoutChoice.available.map(name => (
									<Radio label={name} value={name} key={name} />
								))}
							</RadioGroup>
						</div>
					}
				/>
				{this.renderReembedding()}
			</ButtonGroup>
		);
	}
}

export default Embedding;
